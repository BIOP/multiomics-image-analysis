/**
 * This script is used to detect objects using StarDist in QuPath.
 * It assumes you have fluorescence image.
 *  
 * After defining the builder, it will:
 * 1. Find all selected annotations OR, if no annotations are selected, find all annotations of class TISSUE_CLASS
 * 2. Run the stardist detection using the defined model name or path on the channels defined in CHANNELS_TO_DETECT
 * 3. Create the desired objects with the selected statistics
 * 
 * Please check the documentation and step-by-step tutorial on protocols.io
 * https://go.epfl.ch/multiomics-image-analysis
 * 
 * 
 * NOTE: There are lots of options to customize the detection - this script shows some 
 * of the main ones. Check out other scripts and the QuPath docs for more info.
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

// only one channel
def CHANNELS_TO_DETECT = "DAPI" 

// You can find some at https://github.com/qupath/models
// (Check credit & reuse info before downloading)
def MODEL_PATH = "C:/QuPath_Common_Data_0.7/models/stardist_models/dsb2018_heavy_augment.pb"

def TISSUE_CLASS = "tissue"


/***********************
 * BEGINNING OF THE SCRIPT
 ***********************/


// get the available channels
def imageData = getCurrentImageData()
def originalServer = getCurrentServer()
def originalChannelNames = originalServer.getMetadata().getChannels().collect(e->e.getName())

// get the raw name of the selected channels
def realChannels = []
realChannels.addAll(originalChannelNames.stream().filter(e->e.toLowerCase().contains(CHANNELS_TO_DETECT.toLowerCase())).findAll())

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
println "Starting Stardist..."

// Customize how the StarDist detection should be applied
// Here some reasonable default options are specified
def stardist = StarDist2D
    .builder(MODEL_PATH)
    .channels(realChannels.get(0))            // Extract channel called 'DAPI'
    .normalizePercentiles(1, 99) // Percentile normalization
    .threshold(0.5)              // Probability (detection) threshold
    .pixelSize(0.5)              // Resolution for detection
   // .cellExpansion(5)            // Expand nuclei to approximate cell boundaries
    .measureShape()              // Add shape measurements
    .measureIntensity()          // Add cell measurements (in all compartments)
    .classify(realChannels.get(0)) 
    .build()
	
stardist.detectObjects(imageData, tissueAnnotations)
stardist.close() // This can help clean up & regain memory

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
 
 
import qupath.ext.stardist.StarDist2D
import qupath.lib.scripting.QP