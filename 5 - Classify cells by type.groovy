/**
 * This script classifies cells according to their type, based on their positivity for the defined antibodies.
 * To consider a cell as positive for an antibody, its average intensity has to be larger than the manually-defined
 * threshold for that specific antibody
 * 
 * The cell types are based on the Spyre antibody panels, from Lunaphore
 * https://lunaphore.com/spyre-antibody-panels/. 
 * 
 * Cell classification is split into single and grouped antibodies
 * - single antibodies : the cell is classified as positive for an antibody if its mean intensity is 
 *   larger than the threshold for this antibody
 * 
 * - grouped antibodies : the cells are classified according to their type. If the cell is positive for the
 *   full combinaison of defined antibodies, then it is classified with their respective class.
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


/**********************
 * VARIABLES TO MODIFY
 *********************/
 
 
def antibodyThresholds = [:]
antibodyThresholds["CD45"] = 300
antibodyThresholds["CD3"] = 200
antibodyThresholds["CD4"] = 400
antibodyThresholds["CD8"] = 370
antibodyThresholds["FOXP3"] = 300
antibodyThresholds["CD20"] = 600
antibodyThresholds["CD56"] = 175
antibodyThresholds["CD11c"] = 500
antibodyThresholds["CD68"] = 200
antibodyThresholds["αSMA"] = 750
antibodyThresholds["PD-L1"] = 200
antibodyThresholds["PD-1"] = 400
antibodyThresholds["Ki-67"] = 300
antibodyThresholds["Cytokeratin"] = 1000


def TISSUE_CLASS = "tissue"


/***********************
 * BEGINNING OF THE SCRIPT
 ***********************/
 
  
def singleAntibodies = [
    "αSMA",
    "PD-L1",
    "PD-1",
    "Cytokeratin",
    "Ki-67",
]

def groupedAntibodies = [
    "FOXP3",
    "CD3",
    "CD8",
    "CD4",
    "CD20",
    "CD56",
    "CD11c",
    "CD68",
    "CD45"
]

def correspondanceMap = [:]
correspondanceMap["CD45"] = "Undefined leukocyte"
correspondanceMap["CD3CD4CD45"] = "Helper T Cell"
correspondanceMap["CD3CD8CD45"] = "Cytotoxic T cell"
correspondanceMap["CD3CD45FOXP3"] = "Regulatory T Cell"
correspondanceMap["CD3CD45"] = "Undefined T Cell"
correspondanceMap["CD3CD4CD8CD45"] = "Undefined T Cell"
correspondanceMap["CD3CD4CD45FOXP3"] = "Undefined T Cell"
correspondanceMap["CD3CD8CD45FOXP3"] = "Undefined T Cell"
correspondanceMap["CD3CD4CD8CD45FOXP3"] = "Undefined T Cell"
correspondanceMap["CD11cCD45"] = "Dendritic cell"
correspondanceMap["CD20CD45"] = "B cell"
correspondanceMap["CD45CD56"] = "NK cell"
correspondanceMap["CD45CD68"] = "Macrophage"

def undefinedClass = "Undefined cell"
def unclassifiedClass = "Unclassified cell"


// getting selected annotation or tissue annotations
def selectedObjects = getSelectedObjects()
def tissueAnnotations = []
if(selectedObjects.isEmpty()) {
    tissueAnnotations = getAnnotationObjects().findAll{ it.getPathClass() == getPathClass(TISSUE_CLASS) }
}else{
    tissueAnnotations = selectedObjects.findAll{ it.isAnnotation()}
}

// create new pathClasses if they doesn't exist
def availablePathClasses = QPEx.getQuPath().getAvailablePathClasses();
def newPathClasses = []
newPathClasses.addAll(correspondanceMap.values())
newPathClasses.add(undefinedClass)
newPathClasses.add(unclassifiedClass)

for(def cellClass : newPathClasses){
    def cellPathClass = PathClass.fromString(cellClass);
    if (!availablePathClasses.contains(cellPathClass)) {
        availablePathClasses.add(cellPathClass);
    }
}
QPEx.getQuPath().getProject().setPathClasses(availablePathClasses);

// get cells
def detectionType = "Cell"
def cells  = obj.getChildObjects().findAll { it.isCell() }
if(cells.isEmpty()){
    cells  = obj.getChildObjects().findAll { it.isDetection() && it.getROI() instanceof  PolygonROI}
    detectionType = "Nucleus"
}

// get all measurements
def resultColumns = Results.getAllMeasurements(cells);
def columnNames = resultColumns.getAllNames()

// channel name extraction for processing
def chList = getCurrentServer().getMetadata().getChannels().collect(e->e.getName())

println "Classify cells..."
cells.each{cell->
    def abPositiveList = []

    // check all grouped antibodies channel 
    groupedAntibodies.each{ab->
        def currentChannel = chList.stream().filter(e->e.toLowerCase().startsWith(ab.toLowerCase())).findAll()
        if(!currentChannel.isEmpty()) {
            currentChannel = currentChannel.get(0)
            def metric = detectionType+": "+currentChannel+": Mean"
    
            def abIntensity = resultColumns.getNumericValue(cell, metric)
            if(abIntensity > antibodyThresholds.get(ab)){
                abPositiveList.add(ab)
                cell.measurements.put(ab, 1)
            }else{
                cell.measurements.put(ab, 0)
            }
        } else {
            Logger.warn("The channel '" +ab+ "' is not present in the list of available channels.")
        }
    }

    // classify cells according to positive grouped antibodies
    def abPositive = abPositiveList.sort().join("")
    def cellType = ""
    if(abPositive.isEmpty()){
       cell.setPathClass(PathClass.fromString(unclassifiedClass))
    }else{
        if( correspondanceMap.containsKey(abPositive) ) {
             cell.setPathClass(PathClass.fromString(correspondanceMap.get(abPositive)))
        }else{
             cell.setPathClass(PathClass.fromString(undefinedClass))
        }
    }

    
    // check single antibodies channel 
    singleAntibodies.each{ab->
        def currentChannel = chList.stream().filter(e->e.toLowerCase().startsWith(ab.toLowerCase())).findAll()
        if(!currentChannel.isEmpty()) {
            currentChannel = currentChannel.get(0)
            def metric = detectionType+": "+currentChannel+": Mean"
    
            def abIntensity = resultColumns.getNumericValue(cell, metric)
            if(abIntensity > antibodyThresholds.get(ab)){
                cell.measurements.put(ab, 1)
            }else{
                cell.measurements.put(ab, 0)
            }
        } else {
            Logger.warn("The channel '" +ab+ "' is not present in the list of available channels.")
        }
    }
}

println "End of the script"
return


/***********************
 * IMPORTS
 ***********************/
 
 
import qupath.ext.biop.utils.Results
import qupath.lib.objects.PathObject;
import qupath.lib.gui.scripting.QPEx;