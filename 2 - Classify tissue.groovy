/**
 * This script applies a pixel classifier to detect the tissue
 * Pre-requisit: You need to train your own pixel classifier
 * 
 * Please check the documentation and step-by-step tutorial on protocols.io
 * https://go.epfl.ch/multiomics-image-analysis
 * 
 * 
 * author: Rémy Dornier - PTBIOP
 * date: 2025-10-01
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
 
 
 def CLASSIFIER_NAME = "tissue_classifier"
 def MIN_OBJECT_SIZE = 10000.0
 def MIN_HOLE_SIZE = 5000.0
 def TISSUE_CLASS = "tissue"
 
/***********************
 * BEGINNING OF THE SCRIPT
 ***********************/
 
println "Starting pixel classification..."
createAnnotationsFromPixelClassifier(CLASSIFIER_NAME, MIN_OBJECT_SIZE, MIN_HOLE_SIZE, "SPLIT")

// set a unique name to the tissue annotations
getAnnotationObjects().findAll{ it.getPathClass() == getPathClass(TISSUE_CLASS) }.each {
   it.setName(String.valueOf(it.getID())) 
}

println "End of the script"
return
