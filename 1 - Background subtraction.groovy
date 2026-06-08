/**
 * This script computes the background subtraction on each channel, as it is done in Horizon software.
 * Please check the documentation and step-by-step tutorial on protocols.io
 * 
 * 
 * 
 * author: Rémy Dornier - PTBIOP
 * date: 2025-10-01
 * version: 1.2.0
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
 

// OMERO settings, if your images are stored on OMERO database.
// Otherwise, you don't need to modify them
def HOST = "omero-server.epfl.ch"
def PORT = 4064


/***********************
 * BEGINNING OF THE SCRIPT
 ***********************/
 
 
Date start = new Date()

// Set the ImageServer of the original image
def originalServer = getCurrentServer()
def originalImageData = getCurrentImageData()

// extracting the exposure time of each channel
def exposureTimeMap
if(originalServer instanceof OmeroImageServer) {
    println "Getting exposure times for each channel from OMERO..."
    exposureTimeMap = gettingExposureTimesPerChannelFromOmero(HOST, PORT, originalServer)
}else {
    println "Getting exposure times for each channel from local file..." 
    exposureTimeMap = gettingExposureTimesPerChannel(originalServer)     
}

def rna_auto_fluo_exposure_dict = [:]
def prot_auto_fluo_exposure_dict = [:]
def ch_exposure_dict = [:]
def channels = []
def colors = []
def renamedChannels = []

for (int c = 0; c < originalServer.nChannels(); c++) {
    def currentChannel = originalServer.getChannel(c)
    def currentColor = currentChannel.getColor()
    def currentChannelName = currentChannel.getName()
    def exposureTime = exposureTimeMap.get(c)
    println "Exposure time for channel '"+currentChannelName+"' : "+exposureTime + " s"
    
    // do not apply background subtraction on DAPI channels
    if(currentChannel.getName().toLowerCase().contains("dapi")) {
        channels.add(ColorTransforms.createChannelExtractor(c))
        colors.add(currentColor)
        renamedChannels.add(currentChannelName)
    }else {
        // do not apply background subtraction on auto-fluo channels 
        if(currentChannel.getName().toLowerCase().trim().startsWith("auto")) {
            if(currentChannelName.toLowerCase().contains("rna")) {
                rna_auto_fluo_exposure_dict.put(c, exposureTime)
            }else{  
                prot_auto_fluo_exposure_dict.put(c, exposureTime)
            }
            channels.add(ColorTransforms.createChannelExtractor(c))
            colors.add(currentColor)
            renamedChannels.add(currentChannelName)
        }else {
            // apply background subtraction on all other channels  
            ch_exposure_dict.put(c, exposureTime) 
        }
    }
}

def nRNAGroup = rna_auto_fluo_exposure_dict.size()
def nProtGroup = prot_auto_fluo_exposure_dict.size()

// get the position of the first RNA channel, if exists
def lastAutoRNAPos = originalServer.nChannels() // in case of no RNA, give the lower bounds outside the channels range
def doesContainRNA = false
def rnaSortedPositions = new ArrayList<>(new TreeSet<>(rna_auto_fluo_exposure_dict.keySet()))
   
if(nRNAGroup > 0) {
    lastAutoRNAPos = rnaSortedPositions.get(rnaSortedPositions.size() - 1)
    doesContainRNA = true
}

// get the position of the first Prot auto channel.
def firstAutoProtPos = originalServer.nChannels()
if(nProtGroup > 0) {
    firstAutoProtPos = (new ArrayList(new TreeSet<>(prot_auto_fluo_exposure_dict.keySet()))).get(0)
}

// check the number of RNA cycles to be sure that we get an integer number
def nCycles = 0
if(doesContainRNA) {
    nCycles = (firstAutoProtPos - lastAutoRNAPos - 1) / (2 * nRNAGroup)
    if(nCycles % 1 != 0) {
        Logger.error("RNA cycles - The number of channels "+(firstAutoProtPos - lastAutoRNAPos - 1)+
        " doesn't give an integer number of cycles: "+nCycles+". Stop here the analysis.")
        throw new RuntimeException()
    }
}

println "nAutoRNAchannels: "+nRNAGroup
println "nAutoProtChannels: "+nProtGroup
println "nRNACycles: "+nCycles

println "Compute the background subtracted channels"
for (int c : new ArrayList(new TreeSet<>(ch_exposure_dict.keySet()))) {
    def currentChannel = originalServer.getChannel(c)
    def currentColor = currentChannel.getColor()
    def currentChannelName = currentChannel.getName()
    def currentExposureTime = ch_exposure_dict.get(c)
               
    def autoFluoExposure = 0
    def autoFluoName = ""
    
    // process RNA channels
    if(c > lastAutoRNAPos && c < firstAutoProtPos) {
        def autoRNAChannelPosition = rnaSortedPositions.get((c - lastAutoRNAPos - 1) % nRNAGroup)
        autoFluoExposure = rna_auto_fluo_exposure_dict.get(autoRNAChannelPosition)
        autoFluoName = originalServer.getChannel(autoRNAChannelPosition).getName()
    } else {
     // process prot channels
        for (int key : prot_auto_fluo_exposure_dict.keySet()) {
            autoFluoName = originalServer.getChannel(key).getName()
            def reducedAutoFluoName = autoFluoName.toLowerCase().replace("auto", "").replace("prot", "").trim()
            if(currentChannelName.toLowerCase().contains(reducedAutoFluoName.toLowerCase())){
                autoFluoExposure = prot_auto_fluo_exposure_dict.get(key)
                break
            }
        }
    }
    
    // comupte the correction
    if (autoFluoExposure > 0 && currentExposureTime > 0) {
        def scale = currentExposureTime / autoFluoExposure
        def coefficients = Map.of(
            currentChannelName, 1,
            autoFluoName, -scale,
        )    
        channels.add(ColorTransforms.createLinearCombinationChannelTransform(coefficients))
        renamedChannels.add(currentChannelName + " - BS")
    } else {
        channels.add(ColorTransforms.createChannelExtractor(c)) 
        renamedChannels.add(currentChannelName)
        Logger.error("Background subtraction - Could not find the corresponding auto-fluo channel for channel "+currentChannelName+
        ". Leave the channel without any correction in the final image.")
    }
    colors.add(currentColor)
}
   
println "Create the background subtracted image"
double[] posDummyOffset = new double[1]
posDummyOffset[0] = 1

double[] negDummyOffset = new double[1]
negDummyOffset[0] = -1

// create the background subtracted image server
def newServer = new TransformedServerBuilder(originalServer)
                .applyColorTransforms(channels)
                .normalize(SubtractOffsetAndScaleNormalizer.createWithClipRange(posDummyOffset, null, -Double.MAX_VALUE, Double.MAX_VALUE))
                .normalize(SubtractOffsetAndScaleNormalizer.createWithClipRange(negDummyOffset, null, 0, Double.MAX_VALUE))
                .build() 

def imageData = new ImageData<>(newServer)

Platform.runLater(() -> {
    // display image in the project
    getCurrentViewer().setImageData(imageData)
    
     // set right colors
    QPEx.setChannelColors(imageData, colors.toArray(new Integer[0]));
  
    println "Renaming channels"
    def imageDescription = "Channel coefficients for background subtraction:\n"
    for (int c = 0; c < newServer.nChannels(); c++) {
        def backgroundSubtractedChannelName = newServer.getChannel(c).getName()
        def currentChannelName = renamedChannels.get(c)
        
        // keep track of the coefficients for the background correction in the image description
        imageDescription += currentChannelName + " = " + backgroundSubtractedChannelName + "\n"
    }
    QP.setChannelNames(imageData, renamedChannels.toArray(new String[0]));
    
    // set metadata
    def project = getProject()
    def entry = project.getEntry(imageData)
    entry.getMetadata().put("Type", "Background-subtracted")
    entry.setDescription(imageDescription)
    
    entry = project.getEntry(originalImageData)
    entry.getMetadata().put("Type", "Raw")
})

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
 * METHODS
 ***********************/
 

def gettingExposureTimesPerChannel(server) {
    def exposureTimeMap = new HashMap<>()
    
    // getting file path
    String filePath
    try {
        def uris = server.getURIs()
        filePath = new File(uris[0]).getAbsolutePath()
    } catch(Exception e) {
        filePath = new File(new URI(server.getPath())).getAbsolutePath()
    }
    
    // setting metatada
    def serviceFactory = new ServiceFactory()
    def service = serviceFactory.getInstance(OMEXMLService)
    def meta = service.createOMEXMLMetadata()

    try(def reader = new ImageReader()){
        reader.setMetadataStore(meta)
        reader.setId(filePath) 
         
        // get the first serie
        int series = (server instanceof BioFormatsImageServer) ? server.getSeries() : 0
        reader.setSeries(series)
         
        int nPlanes = meta.getPlaneCount(series)
         
        // get the list of exposure times per channel
        if (nPlanes == 0) {
            Logger.warn("No metadata available for this image. No exposure time found")
        } else {
        
            def seen = [] as Set
        
            for (int p = 0; p < nPlanes; p++) {
                def expTime = meta.getPlaneExposureTime(series, p)
                if (expTime != null) {
                    int c = meta.getPlaneTheC(series, p)?.numberValue?.intValue() ?: p
                    if (!seen.contains(c)) {
                        seen.add(c)
                        def chName = server.getMetadata().getChannels()[c].getName()
                        def sec = expTime.value(UNITS.SECOND)
                        exposureTimeMap.put(c, sec)
                    }
                }
            }
        
            if (seen.isEmpty()) {
                Logger.warn("No exposure times found in the metadata.")
            }
        }
    }catch(Exception e) {
       Logger.error("Error during opening the file" +e) 
    }
    
    return exposureTimeMap
}


def gettingExposureTimesPerChannelFromOmero(HOST, PORT, imageServer) {
    def gateway = new Gateway(new IceLogger())
    def ctx = null
    def imageId = imageServer.getId() 
    Unit<ome.units.quantity.Time> unit = SECOND
    def exposureTimeMap = new HashMap<>()

     // connect to OMERO via session ID
    def sessionId = imageServer.getClient().getApisHandler().getSessionUuid().get()
    def cred = new LoginCredentials(sessionId, sessionId, HOST, PORT)
    def connectedUser = gateway.connect(cred);
           
    if(gateway.isConnected()) {
       try{
            println "Connected to "+HOST
            
            // set the security context
            ctx = new SecurityContext(connectedUser.getGroupId())
            ctx.setExperimenter(connectedUser);
            ctx.setServerInformation(cred.getServer());
            
            // get the current image from OMERO
            def omeroImageData = gateway.getFacility(BrowseFacility.class).getImage(ctx, imageId)
            
            // list exposure time per channel
            def pixelData = omeroImageData.getDefaultPixels()
            def planeList = gateway.getFacility(MetadataFacility.class).getPlaneInfos(ctx, pixelData)
            for (def p : planeList) {
               ome.units.quantity.Time t = convertTime(p.getExposureTime());
               if (t != null) {
                   Number value = t.value(unit);
                   if (value != null) {
                       def exposureTime = 0
                       if(exposureTimeMap.containsKey(p.getTheC())) {
                           exposureTime = exposureTimeMap.get(p.getTheC())
                           exposureTime = (exposureTime + value.doubleValue()) / 2
                       }else {
                           exposureTime = value.doubleValue()
                       }
                       exposureTimeMap.put(p.getTheC(), exposureTime)
                   }
               }
            }
        } finally {
           gateway.disconnect()
           println "Disconnected from "+HOST
           if(ctx != null) {
               ctx = new SecurityContext(-1)
               ctx.setExperimenter(new ExperimenterData());
           }
       }
   }
   return exposureTimeMap
}


/***********************
 * IMPORTS
 ***********************/
 
 
import java.util.stream.IntStream
import java.util.stream.Stream
import qupath.lib.images.servers.*
import qupath.ext.omero.core.imageserver.OmeroImageServer
import qupath.ext.omero.core.pixelapis.ice.IceLogger
import omero.gateway.Gateway;
import omero.gateway.LoginCredentials;
import omero.gateway.SecurityContext;
import omero.gateway.model.ExperimenterData;
import omero.gateway.facility.BrowseFacility;
import omero.gateway.facility.MetadataFacility;
import omero.model.Time;
import omero.model.TimeI;
import ome.formats.model.UnitsFactory;
import ome.units.unit.Unit;
import static java.lang.Double.NaN;
import static ome.formats.model.UnitsFactory.convertTime;
import static ome.units.UNITS.SECOND;
import ome.units.UNITS
import java.awt.Color;
import qupath.lib.gui.scripting.QPEx;
import qupath.lib.images.servers.transforms.SubtractOffsetAndScaleNormalizer
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import qupath.lib.images.servers.bioformats.BioFormatsImageServer
import loci.formats.ImageReader
import loci.common.services.ServiceFactory
import loci.formats.services.OMEXMLService