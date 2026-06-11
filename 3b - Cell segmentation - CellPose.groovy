/**
 * This script is used to detect objects using a Cellpose model within QuPath.
 * After defining the builder, it will:
 * 1. Find all selected annotations OR, if no annotations are selected, find all annotations of class TISSUE_CLASS
 * 2. Export the annotations to a temp folder that can be specified with tempDirectory()
 * 3. Run the cellpose detection using the defined model name or path on the channels defined in CHANNELS_TO_DETECT
 * 4. Reimport the mask images into QuPath and create the desired objects with the selected statistics
 *
 * Please check the documentation and step-by-step tutorial on protocols.io
 * https://go.epfl.ch/multiomics-image-analysis
 * 
 * 
 * NOTE: that this template does not contain all options, but should help get you started
 * See all options in https://biop.github.io/qupath-extension-cellpose/qupath/ext/biop/cellpose/CellposeBuilder.html
 * and in https://cellpose.readthedocs.io/en/latest/command.html
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


def CHANNELS_TO_DETECT = [
    "DAPI",
    "CD3", 
    "FoxP3"
]

// Specify the model name (cyto, nuclei, cyto2, ... or a path to your custom model as a string)
// Other models for Cellpose https://cellpose.readthedocs.io/en/latest/models.html
def MODEL_PATH = "cpsam"

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

Date start = new Date()
println "Starting CellPose..."

def cellpose = Cellpose2D.builder( MODEL_PATH )
        .pixelSize( 0.5 )                      // Resolution for detection in um
        .channels( realChannels as String[] )	               // Select detection channel(s)
//        .tempDirectory( new File( '/tmp' ) )         // Temporary directory to export images to. defaults to 'cellpose-temp' inside the QuPath Project
//        .preprocess( ImageOps.Filters.median( 1 ) )  // List of preprocessing ImageOps to run on the images before exporting them
        .normalizePercentilesGlobal( 0.1, 99.8, 10 ) // Convenience global percentile normalization. arguments are percentileMin, percentileMax, dowsample. If no channel mentioned, arguments are applied for all selected channels. Otherwise, applied only for the selected channel. See https://forum.image.sc/t/cellpose-sam-qupath-extension/112418/27
//        .tileSize( 1024 )                  // If your GPU can take it, make larger tiles to process fewer of them. Useful for Omnipose
//        .cellposeChannels( 1,2 )           // Overwrites the logic of this plugin with these two values. These will be sent directly to --chan and --chan2
//        .cellprobThreshold( 0.0 )          // Threshold for the mask detection, defaults to 0.0
        .flowThreshold( 1 )              // Threshold for the flows, defaults to 0.4
//        .diameter( 15 )                    // Median object diameter. Set to 0.0 for the `bact_omni` model or for automatic computation
//        .useOmnipose()                     // Use omnipose instead
        .useCellposeSAM()                  // Use cellposeSAM (i.e. cellpose 4.x.x) env instead of previous versions of cellpose <= v3.x.x
//        .addParameter( "cluster" )         // Any parameter from cellpose or omnipose not available in the builder.
//        .addParameter( "save_flows" )      // Any parameter from cellpose or omnipose not available in the builder.
//        .addParameter( "anisotropy", "3" ) // Any parameter from cellpose or omnipose not available in the builder.
        .cellExpansion( 3 )              // Approximate cells based upon nucleus expansion
//        .cellConstrainScale( 1.5 )         // Constrain cell expansion using nucleus size
//        .classify( "My Detections" )       // PathClass to give newly created objects
        .measureShape()                    // Add shape measurements
        .measureIntensity()                // Add cell measurements (in all compartments)
//        .createAnnotations()               // Make annotations instead of detections. This ignores cellExpansion
//        .simplify( 0 )                     // Simplification 1.6 by default, set to 0 to get the cellpose masks as precisely as possible
//        .disableGPU()                      // Force using CPU.
//        .excludeEdges()                    // remove cells touching the border => higher priority than constrainToParent()
//        .constrainToParent(false, 15)      // display all and entirely the cells intersecting the parent annotation, with an optional padding around the parent annotation given in um. Default true and 15um. Ignored if excludeEdges() is called.
        .build()


cellpose.detectObjects( imageData, tissueAnnotations )

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


import qupath.ext.biop.cellpose.Cellpose2D