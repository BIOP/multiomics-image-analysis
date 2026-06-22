/**
 * This script is used to detect objects using a Spotiflow model within QuPath.
 * After defining the builder, it will:
 * 1. Find all selected annotations OR, if no annotations are selected, find all annotations of class TISSUE_CLASS
 * 2. Export the selected annotations to a temp folder that can be specified with tempDirectory()
 * 3. Run the spotiflow detction using the defined/default model name or path on the channels defined in CHANNELS_TO_DETECT
 * 4. Create the desired objects (i.e. points) with the selected statistics (i.e. spotiflow outputs)
 *
 * NOTE: that this template does not contain all options, but should help get you started
 * See all options by calling spotiflow.helpPredict()
 *
 * Please check the documentation and step-by-step tutorial on protocols.io
 * https://go.epfl.ch/multiomics-image-analysis
 * 
 * 
 * author: Rémy Dornier - PTBIOP
 * date: 2026-06.08
 * version: 1.0.0
 * 
 * Last tested on QuPath 0.7.x
 * 
 * -----------------------------------------------------------------------------
 * Copyright (c) 2026 ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, BioImaging And Optics Platform (BIOP)
 * All rights reserved.
 * 
 * Licensed under the BSD-3-Clause License:
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided 
 * that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
 *    in the documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products 
 *     derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
 * BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * -----------------------------------------------------------------------------
 */


/***********************
 * VARIABLES TO MODIFY
 ***********************/

// If you have trained a custom model, specify the model directory as a File in setModelDir()
// If you want to use any other pre-trained models, specify its name in setPretrainedModelName()
// -> List of all pre-trained models : https://weigertlab.github.io/spotiflow/pretrained.html
def PRETRAINED_MODEL = "general"

// no channel limit
def CHANNELS_TO_DETECT = [
    "T5",
    "T6",
    "T7",
]

def TISSUE_CLASS = "tissue"


/***********************
 * BEGINNING OF THE SCRIPT
 ***********************/


// get the available channels
def originalServer = getCurrentServer()
def originalChannelNames = originalServer.getMetadata().getChannels().collect(e->e.getName())
def imageData = getCurrentImageData()

// get the raw name of the selected channels
def realChannels = []
CHANNELS_TO_DETECT.each {
   realChannels.addAll(originalChannelNames.stream().filter(e->e.toLowerCase().contains(it.toLowerCase())).findAll())
}

// remove duplicate channels
realChannels = realChannels.unique()

// check if selected channels are found
if(realChannels.isEmpty()) {
   Logger.warn("The channels you give cannot be retrieved in the list of available channels. Please check the channel names") 
   Logger.warn("Channels: "+ CHANNELS_TO_DETECT) 
   return
}else {
   Logger.info("The following channels will be used for detection:")
   Logger.info("Channels: "+ realChannels) 
}


// get tissue annotations
def selectedObjects = getSelectedObjects()
def tissueAnnotations = []
if(selectedObjects.isEmpty()){
    tissueAnnotations = getAnnotationObjects().findAll{ it.getPathClass() == getPathClass(TISSUE_CLASS) }
}else{
    tissueAnnotations = selectedObjects.findAll{ it.isAnnotation()}
}
println("Found "+tissueAnnotations.size()+" "+TISSUE_CLASS+" annotation(s)")

if (tissueAnnotations.isEmpty()) {
    Logger.error("No parent objects are selected!")
    return
}

// set a unique name to the tissue annotation if it's not already the case
tissueAnnotations.each {
    if(it.getName() == null || it.getName().isEmpty()) {
        it.setName(String.valueOf(it.getID()))  
    }
}   

Date start = new Date()
println "Starting Spotiflow..."

def spotiflow = Spotiflow.builder()
//        .tempDirectory(new File("path/to/tmp/folder"))       // OPTIONAL : default is in 'qpProject/spotiflow-temp' folder
//        .setModelDir(new File("path/to/my/model"))           // OPTIONAL : path to your own trained model
        .setPretrainedModelName(PRETRAINED_MODEL)                 // OPTIONAL : Default is 'general'
//        .setMinDistance(2)                                   // OPTIONAL : Positive integer value
        .setProbabilityThreshold(0.6)                        // OPTIONAL : Positive value
//        .disableGPU()                                        // OPTIONAL : Force using CPU ; default is automatic (let spotiflow decide)
//        .process3d()                                         // OPTIONAL : process the entire zstack
//        .zPositions(0,5)                                     // OPTIONAL : ONLY works wih process3d(). Select a sub-stack (start and end inclusive)
//        .doSubpixel(true)                                    // OPTIONAL : true to get subpixel resolution ; false to not. Default: let spotiflow choose
//        .setClass("ClassName")                               // OPTIONAL : set the same class for all detections. Default: not assign any classes
        .setClassChannelName()                               // OPTIONAL : create a new class for each channel and assign detection to it. Default: not assign any classes
//        .nThreads(12)                                        // OPTIONAL : How much you want to paralellize processing. Default 12
//        .saveBuilder("MyFancyName")                          // OPTIONAL : To save builder parameters as JSON file
        .saveTempImagesAsOmeZarr()                           // OPTIONAL : ONLY AVAILABLE FOR SPOTIFLOW >= 0.5.8. Save temp images as ome-zarr instead of ome.tiff
//        .clearAllChildObjects()                              // OPTIONAL : Clear all previous detections, whatever their class
//        .createAnnotations()                                 // OPTIONAL : Create annotations instead of detections. WARNING: this can slow up a lot QuPath. Only to use to pre-annotated small patches for later training.
        .clearChildObjectsBelongingToCurrentChannels()       // OPTIONAL : Clear all previous detections which belong to the current selected channels (i.e. with their class set with the name of the channel)
        .channels(realChannels as String[])        // REQUIRED : list of channel name(s) to process. At least one channel is required
        .cleanTempDir()                                      // OPTIONAL : Clean all files from the tempDirectory
//        .addParameter("key","value")                         // OPTIONAL : Add more parameter, base on the available ones
        .build()


// print the available arguments for prediction
//spotiflow.helpPredict()

// detect spots
spotiflow.detectObjects( imageData, getProjectEntry().getID(), tissueAnnotations )

// prints processing time
Date stop = new Date()
long milliseconds = stop.getTime() - start.getTime()
int seconds = (int) (milliseconds / 1000) % 60 ;
int minutes = (int) ((milliseconds / (1000*60)) % 60);
int hours   = (int) ((milliseconds / (1000*60*60)) % 24);
println "Processing done in " + hours + " hour(s) " + minutes + " minute(s) " + seconds + " second(s)"

println "End of the script"
return


/***********************
 * IMPORTS
 ***********************/
 
 
import qupath.ext.biop.spotiflow.Spotiflow