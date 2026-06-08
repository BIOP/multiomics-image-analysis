/**
 * This script is used to detect objects using a InstanSed model within QuPath.
 * After defining the builder, it will:
 * 1. Find all selected annotations OR, if no annotations are selected, find all annotations of class TISSUE_CLASS
 * 2. Run the instanseg detection using the defined model name or path on the channels defined in CHANNELS_TO_DETECT
 * 3. Create the desired objects with the selected statistics
 * 
 * Before running the script, please go under "Extensions > Deep Java Library > Manage DJL engines"
 * and click on PyTorch > check/download
 * 
 * NOTE: If the instenseg model is located on a server (i.e. \\sv-nas1.rcp....), you'll receive an error.
 * To get rid of this, save your project, close it, move it to your local computer et start it again.
 *  
 * author: Rémy Dornier - PTBIOP
 * date: 2026-01-19
 * version: 1.0.0
 * 
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


/**********************
 * VARIABLES TO MODIFY
 *********************/


def DEVICE = "gpu" // change to cpu if no gpu available
def MODEL_PATH = "C:/QuPath_Common_Data_0.7/models/instanseg_models/fluorescence_nuclei_and_cells-0.1.0"

def CHANNELS_TO_DETECT = [
    "DAPI",
    "FITC", 
    "Cy3", 
    "Cy5", 
]

def TISSUE_CLASS = "tissue"

/***********************
 * BEGINNING OF THE SCRIPT
 ***********************/


// instanseg parameters
def threads = 10
def tileSize = 1024
def padding = 64
def makeMeasurements = true
def randomColors = false


// get the annotations to process
println "Getting objects to process"

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

// channel name extraction for processing
def chList = getCurrentServer().getMetadata().getChannels().collect(e->e.getName())
restrictedChList = []
CHANNELS_TO_DETECT.each {
   restrictedChList.addAll(chList.stream().filter(e->e.toLowerCase().contains(it.toLowerCase())).findAll())
}

//Create the colorTransforms for instanseg.
def instansegColorTransforms = []
restrictedChList.each{
    instansegColorTransforms.add(ColorTransforms.createChannelExtractor(it))
}

// select the annotations to process
selectObjects(tissueAnnotations)

// run instanseg
println "Running instanseg on channels : " + restrictedChList
Date start = new Date()

qupath.ext.instanseg.core.InstanSeg.builder()
    .modelPath(MODEL_PATH)
    .device(DEVICE)
    .nThreads(threads)
    .tileDims(tileSize)
    .interTilePadding(padding)
    .inputChannels(instansegColorTransforms)
    .outputChannels()
    .makeMeasurements(makeMeasurements)
    .randomColors(randomColors)
    .build()
    .detectObjects()


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
 
 
import qupath.lib.images.servers.ColorTransforms
import qupath.lib.gui.scripting.QPEx  
import java.util.stream.Collectors