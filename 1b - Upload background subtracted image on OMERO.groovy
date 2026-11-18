 /**
 * This script aims at importing the background subtracted image, generated from the 
 * multiomics image analysis protocol in QuPath, on OMERO.
 *
 * Please check the documentation and step-by-step tutorial on protocols.io
 * https://go.epfl.ch/multiomics-image-analysis
 * 
 * 
 * dependencies
 *     - raw2ometiff
 *     - qupath-extension-omero
 * 
 * 
 * author: Rémy Dornier - PTBIOP
 * date: 2026.11.16
 * version: 1.3.0
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
 
 
// OMERO parameters
def HOST = "omero-server.epfl.ch"
def PORT = 4064

// OME-Zarr export parameters
def N_THREADS = 64
def TILE_SIZE = 1024

// Path to the NGFF converter
// For Windows, file should be located under /path/to/raw2ometiff-x.x.x/bin/raw2ometiff.bat
// For Mac/Linux, file should be located under /path/to/raw2ometiff-x.x.x/bin/raw2ometiff
def RAW2OMETIFF_EXE = "C:/bfconvert/raw2ometiff-0.10.1/bin/raw2ometiff.bat"

// Temporary output folder where to save the ome.zarr and the ome.tiff
def OUTPUT_FOLDER = "path/to/results_folder"

// Boolean to delete intermediate results
// true = delete
// false = don't delete
def DELETE_OME_ZARR_TEMP = true
def DELETE_OME_TIFF_TEMP = false


/***********************
 * BEGINNING OF THE SCRIPT
 ***********************/


// get current image properties
def imageEntry = getProjectEntry()
def imageName = imageEntry.getImageName()
def imageData = imageEntry.readImageData()
def imageServer = imageData.getServer()

// create the full image annotation
def parentRoi = ROIs.createRectangleROI(0, 0, imageServer.getWidth(), imageServer.getHeight())
def region = RegionRequest.createInstance(imageData.getServerPath(), 1.0, parentRoi);

// make sure that the current image is the background subracted one
def currentImageType = getProjectEntry().getMetadata().get("Type")
if(!currentImageType.equalsIgnoreCase("background-subtracted")) {
   Logger.error("The current image is not a background subtracted image. Please open a background subtracted image and re-run the script") 
   return
}

// make sure that the raw image is part of the current QuPath project
def rawEntry = null
for(def entry: getProject().getImageList()){
    if(entry.getMetadata().get("Type") != null && entry.getMetadata().get("Type").equalsIgnoreCase("raw") && imageName.contains(entry.getImageName())) {
        rawEntry = entry
        break
    }
}

if(rawEntry == null) {
   Looger.error("The corresponding raw image cannot be found in the QuPath project. Cannot save the image back to OMERO")
   return
}
    
// convert the background subracted image to OME-Zarr
println "Exporting background subtracted image to OME.Zarr..."
def omeZarrPath = saveOmeZarrFile(imageServer, N_THREADS, TILE_SIZE, region, OUTPUT_FOLDER, imageName)

// make sure that the 'Big Endian' property is part of the Metadata.ome.xml
println "Reparing OME.Zarr..."
repareOmeZarrForOmeroCompatibility(omeZarrPath)

// build the ome.tiff path
def omeZarrParentPath = new File(omeZarrPath).getParent()
def outputTiff = omeZarrParentPath + File.separator + imageName + ".ome.tiff"

// convert the ome.zarr into ome.tiff file
println "Converting background subtracted image to OME.Tiff file..."
convertZarrToOmeTiff(RAW2OMETIFF_EXE, omeZarrPath, outputTiff, N_THREADS)

// import background subtracted ome.tiff image on OMERO
println "Importing background subtracted image on OMERO..."
def ids = importToOmero(HOST, PORT, rawEntry, outputTiff)
if(!ids.isEmpty()) {
    deleteTempFiles(DELETE_OME_ZARR_TEMP, DELETE_OME_TIFF_TEMP, omeZarrPath, outputTiff) 
    
    println "Sending channel colors & display settings to OMERO"
    sendingImageSettings(imageData, imageServer, rawEntry, ids.get(0))
}

println "End of script"
return


/***********************
 * METHODS
 ***********************/
 
 
def saveOmeZarrFile(imageServer, nThreads, tileSize, region, outputFolder, imageName) {
    
    // get the list of downsample levels
    List<Double> downsampleList = new ArrayList<>();
    double nextDownsample = 1.0;
    do {
        downsampleList.add(nextDownsample);
        nextDownsample *= 2;
    }while ((int)(imageServer.getWidth() / nextDownsample) > tileSize && (int)(imageServer.getHeight() / nextDownsample) > tileSize);
    double[] downsamples = downsampleList.stream().mapToDouble(d -> d).toArray();


    // building ome.zarr writer
    def outputPath = outputFolder + File.separator + imageName + ".ome.zarr"
    logger.info("Saving image(s) into the temporary folder as OME-Zarr : {}", outputPath);
    OMEZarrWriter.Builder builder = new OMEZarrWriter.Builder(imageServer);
    builder.parallelize(nThreads)
            .tileSize(tileSize)
            .region(region)
            .downsamples(downsamples);
    
    // save ome-zarr
    try {
        OMEZarrWriter omeZarrWriter = builder.build(outputPath);
        omeZarrWriter.writeImage();
    } catch (Exception e) {
        logger.error("Error during writing OME-Zarr file", e);
    }
    return outputPath
}


def repareOmeZarrForOmeroCompatibility(omeZarrPath) {
    
    // get the path to the Zarr Metadata file
    def xmlFile = new File(omeZarrPath + File.separator + "OME" + File.separator + "METADATA.ome.xml")

    // rename the original file to keep it clean in case there is an issue during modification
    def xmlFileRename = new File(xmlFile.getAbsolutePath() + ".beforeModification")
    FileUtils.copyFile(xmlFile, xmlFileRename);
    
    // Create a DocumentBuilder
    DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
    DocumentBuilder builder = factory.newDocumentBuilder();
    
    // Parse the XML file
    Document document = builder.parse(xmlFile);
    
    // Access Pixels by tag name
    NodeList pixelsList = document.getElementsByTagName("Pixels")
    Node pixels = pixelsList.item(0)
    if (pixels.getNodeType() == Node.ELEMENT_NODE) {
        Element pixelElement = pixels
        
        // add the tag BigEndian if it doesn't exist
        if(pixelElement.getAttribute("BigEndian").isEmpty()){
        	pixelElement.setAttribute("BigEndian", "false");
        }
    }
    
    // Produce string XML
    try(OutputStream outputStream = new FileOutputStream(xmlFile)){
        outputStream.write(asString(document).getBytes());
    } catch(Exception e){
        e.printStackTrace();
        throw e;
    }
}


def convertZarrToOmeTiff(raw2ometiffExe, omeZarrPath, outputTiff, nThreads) {
    // define command line arguments
    def extraArgs = ["--overwrite", "--compression", "LZW", "--max_workers", String.valueOf(nThreads)] // too much info logs, but in case: "--log-level", "INFO"
    
    List<String> shell = []
    
    // choose the right OS command
    switch (org.controlsfx.tools.Platform.getCurrent()) {
        case org.controlsfx.tools.Platform.UNIX:
        case org.controlsfx.tools.Platform.OSX:
            shell.addAll(Arrays.asList("bash", "-c"));

            // If there are spaces, then we should encapsulate the command with quotes
            command = command.stream().map(s -> {
                        if (s.trim().contains(" "))
                            return "\"" + s.trim() + "\"";
                        return s;
                    }).collect(Collectors.toList());

            // The last part needs to be sent as a single string, otherwise it does not run
            String cmdString = command.toString().replace(",","");

            shell.add(cmdString.substring(1, cmdString.length()-1));
            break;

        case org.controlsfx.tools.Platform.WINDOWS:
        default:
            shell.addAll(Arrays.asList("cmd.exe", "/C"));
            break;
    }
    
    // build the shell command
    shell << raw2ometiffExe
    shell << omeZarrPath
    shell << outputTiff
    shell.addAll(extraArgs)

    println "Full command : " + shell.join(" ")

    try {
    	// start the ome.tiff conversion with the command line
        ProcessBuilder pb = new ProcessBuilder(shell)
        pb.redirectErrorStream(true)
        Process process = pb.start()
        
        // Keep the log of the process
        def logResults = new ArrayList<>();
        Thread t = new Thread(Thread.currentThread().getName() + "-" + process.hashCode()) {
            @Override
            public void run() {
                BufferedReader stdIn = new BufferedReader(new InputStreamReader(process.getInputStream()));
                try {
                    for (String line = stdIn.readLine(); line != null; ) {
                        Logger.info("{}: {}", name, line);
                        logResults.add(line);
                        line = stdIn.readLine();
                    }
                } catch (IOException e) {
                    Logger.warn(e.getMessage());
                }
            }
        };
        t.setDaemon(true);
        t.start();
        
        // wait for the end of the ome.tiff conversion
        Logger.info("Virtual Environment Runner Started");
        int exitCode = process.waitFor()
    
    	// log output
        if (exitCode == 0) {
            Logger.info("Successful conversion from ome.zarr to ome-tiff : ${outputTiff}")
        } else {
            Logger.error("Issue during the process (return code : ${exitCode})")
        }
    
    } catch (Exception e) {
        Logger.error("Exception during raw2ometiff execution : ${e.getMessage()}")
        e.printStackTrace()
    }   
}


def importToOmero(host, port, rawEntry, omeTiffPath) {
    // read the current QuPath OMERO session
    def sessionId = rawEntry.readImageData().getServer().getClient().getApisHandler().getSessionUuid().get()
    def gateway = new Gateway(new IceLogger())
    def cred = new LoginCredentials(sessionId, sessionId, host, port)
    def connectedUser = gateway.connect(cred);
    def ctx = null
    def success = false 
    
    if(gateway.isConnected()) {
       try{
           Logger.info("Connected to "+host)
           ctx = new SecurityContext(connectedUser.getGroupId())
           ctx.setExperimenter(connectedUser);
           ctx.setServerInformation(cred.getServer());
           
           // get the target dataset id
           def datasetId = Long.parseLong(rawEntry.getMetadata().get("dataset-id"))
           def datasets = gateway.getFacility(BrowseFacility.class).getDatasets(ctx, Collections.singletonList(datasetId))
           if (datasets.isEmpty()) {
                Logger.error("Dataset " + id+ " doesn't exist in this context");
                throw new NoSuchElementException(msg);
           }
           def dataset = datasets.get(0)
           
           // import background subtracted image on OMERO
           success = importImages(gateway, ctx, dataset, 1, omeTiffPath)
           if(success) {
              Logger.info("Successful import of background image on OMERO.") 
           }else {
              Logger.error("Background subtracted image wasn't imported on OMERO ; an issue occurred duing the import.") 
           }
           
       }finally {
       	// disconnect from OMERO
           gateway.disconnect()
           Logger.info("Disconnected from "+host)
           if(ctx != null) {
               ctx = new SecurityContext(-1)
               ctx.setExperimenter(new ExperimenterData());
           }
       }
    }
    return success
}


def importImages(client, ctx, target, threads, path){
    // setup the import
    ImportConfig config = new ImportConfig();
    String type = PojoMapper.getGraphType(target.getClass());
    def user = ctx.getExperimenterData()
    config.username.set(user.getUserName());
    config.email.set(user.getEmail());
    Collection<Pixels> pixels = new ArrayList<>(1);

	// open OMERO import connection
    OMEROMetadataStoreClient store = client.getImportStore(ctx);
    try (OMEROWrapper reader = new OMEROWrapper(config)) {
        store.logVersionInfo(config.getIniVersionNumber());
        reader.setMetadataOptions(new DefaultMetadataOptions(MetadataLevel.ALL));

        ImportLibrary library = new ImportLibrary(store, reader);
        library.addObserver(new LoggingImportMonitor());

        ErrorHandler handler = new ErrorHandler(config);

        ImportCandidates candidates = new ImportCandidates(reader, new String[]{path}, handler);

        ExecutorService uploadThreadPool = Executors.newFixedThreadPool(config.parallelUpload.get());
		
		// import !
        List<ImportContainer> containers = candidates.getContainers();
        if (containers != null) {
            for (int i = 0; i < containers.size(); i++) {
                ImportContainer container = containers.get(i);
                container.setTarget(target.asIObject());
                List<Pixels> imported = library.importImage(container, uploadThreadPool, i);
                pixels.addAll(imported);
            }
        }
        uploadThreadPool.shutdown();
    } catch (Throwable e) {
        throw new RuntimeException(e);
    } finally {
        client.closeImport(ctx, null);
    }
	
	// return the imported ids
    List<Long> ids = new ArrayList<>(pixels.size());
    pixels.forEach(pix -> ids.add(pix.getImage().getId().getValue()));
    return ids.stream().distinct().collect(Collectors.toList());
}


def deleteTempFiles(deleteOmeZarrTemp, deleteOmeTiffTemp, omeZarrPath, outputTiff) {
   // delete the ome.zarr file
   if(deleteOmeZarrTemp) {
   	  println "Deleting temporary OME.Zarr.."
      def omeZarrFile = new File(omeZarrPath) 
      if(omeZarrFile.exists()){
          FileUtils.deleteDirectory(omeZarrFile);
      }
   }
   
   // delete the ome.tiff file
   if(deleteOmeTiffTemp) {
   	  println "Deleting temporary OME.Tiff.."
      def omeTiffFile = new File(outputTiff) 
      if(omeTiffFile.exists()){
          omeTiffFile.delete()
      }
   }   
}


def sendingImageSettings(imageData, imageServer, rawEntry, imageId) {
    def apiHandler = rawEntry.readImageData().getServer().getClient().getApisHandler()
    
    // Retrieve channel colors from current image
    def channelColors = imageServer.getMetadata().getChannels().stream().map(ImageChannel::getColor).toList()

    // Send channel colors
    apiHandler.changeChannelColors(imageId, channelColors).get()

    println "Channel colors sent"
   

    // Retrieve display ranges from current image
    def displayRanges = getCurrentViewer().getImageDisplay().availableChannels().stream()
            .map(channel -> new ChannelSettings(channel.getMinDisplay(), channel.getMaxDisplay()))
            .toList()

    // Send channel display ranges
    apiHandler.changeChannelDisplayRanges(imageId, displayRanges).get()

    println "Channel display ranges sent"
}


def asString(Document document) {
    TransformerFactory transformerFactory = TransformerFactory.newInstance();
    Transformer transformer = transformerFactory.newTransformer();

    //Setup indenting to "pretty print"
    transformer.setOutputProperty(OutputKeys.INDENT, "yes");
    transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "4");

    Source source = new DOMSource(document);
    ByteArrayOutputStream os = new ByteArrayOutputStream();
    Result result = new StreamResult(new OutputStreamWriter(os, "utf-8"));
    transformer.transform(source, result);

    return os.toString();
}


/***********************
 * IMPORTS
 ***********************/
 
 
import org.apache.commons.io.FileUtils;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import org.w3c.dom.Document;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import javax.xml.transform.OutputKeys;
import org.w3c.dom.NodeList;
import javax.xml.transform.Source;
import org.w3c.dom.Element;
import javax.xml.transform.Result;
import org.w3c.dom.Node;
import java.io.File;
import qupath.lib.images.writers.ome.zarr.OMEZarrWriter;
import qupath.ext.omero.core.pixelapis.ice.IceLogger
import omero.gateway.Gateway;
import omero.gateway.LoginCredentials;
import omero.gateway.SecurityContext;
import omero.gateway.model.ExperimenterData;
import loci.formats.in.DefaultMetadataOptions;
import loci.formats.in.MetadataLevel;
import ome.formats.OMEROMetadataStoreClient;
import ome.formats.importer.ImportCandidates;
import ome.formats.importer.ImportConfig;
import ome.formats.importer.ImportContainer;
import ome.formats.importer.ImportLibrary;
import ome.formats.importer.OMEROWrapper;
import ome.formats.importer.cli.ErrorHandler;
import ome.formats.importer.cli.LoggingImportMonitor;
import omero.ServerError;
import omero.gateway.model.DataObject;
import omero.gateway.util.PojoMapper;
import omero.model.Pixels;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import omero.gateway.facility.BrowseFacility;
import java.util.stream.Collectors
import qupath.ext.omero.core.apis.commonentities.ChannelSettings
import qupath.ext.omero.core.imageserver.OmeroImageServer
import qupath.lib.images.servers.ImageChannel