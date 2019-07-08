"#S3segmenter" 
S3segmenter is a Matlab-based set of functions that generates single cell (nuclei and cytoplasm) label masks. Inputs are:
1. an .ome.tif (preferably flat field corrected)
2. a 3-class probability maps derived from a deep learning model such as UNet. Classes include background, nuclei contours, and nuclei foreground.

The centers of each nuclei are obtained by finding local maxima from the nuclei foreground. These are used for marker-controlled watershed constrained by the nuclei contours. 

To segment cytoplasm, the nuclei are in turn used for a marker-controlled watershed segmentation constrained by a cytoplasmic marker such as B-catenin. The channel number of this marker must be specified. A 3-pixel annulus around each nucleus will also be used to segment cytoplasm.

How to run:
In Matlab, set path to the folder of the cloned repo. Type:
O2batchS3segmenterWrapperR('/path/to/files/')

Use the following name-value pairs arguments to customize the code to your experiment:
ip.addParamValue('HPC','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('fileNum',1,@(x)(numel(x) > 0 & all(x > 0 )));  % if using a cluster, this specifies which file index to work on
ip.addParamValue('CytoMaskChan',[2],@(x)(numel(x) > 0 & all(x > 0 )));  % select any number of channels for cytoplasm
ip.addParamValue('TissueMaskChan',[3],@(x)(numel(x) > 0 & all(x > 0 ))); % select any number of channels for tissue mask
ip.addParamValue('RefineTissueMask',[0],@(x)(numel(x) > 0 & all(x > 0 ))); % constrict the tissue mask to eliminate high autofluorescent regions

ip.addParamValue('mask','tissue',@(x)(ismember(x,{'TMA','tissue','none'}))); % set to true if sample is TMA cores
ip.addParamValue('crop','noCrop',@(x)(ismember(x,{'interactiveCrop','autoCrop','dearray','noCrop'}))); % interactiveCrop - a GUI-based crop selector, 'autoCrop' - takes the middle third region,'dearray', set to true if using TMA cores, 'noCrop', no cropping

ip.addParamValue('cytoMethod','distanceTransform',@(x)(ismember(x,{'RF','distanceTransform','bwdistanceTransform','ring'})));
ip.addParamValue('nucleiFilter','IntPM',@(x)(ismember(x,{'LoG','Int','IntPM','none'}))); % feature to threshold nuclei. 'IntPM' - intensity of probability map, 'Int' - intensity of DAPI channel, 'LoG', intensity of LoG filter response, 'none', accept all nuclei

ip.addParamValue('measureFeatures','false',@(x)(ismember(x,{'true','false'}))); % extracts intensity features from mask
ip.addParamValue('nucleiRegion','watershedContourInt',@(x)(ismember(x,{'watershedContourDist','watershedContourInt','watershedBWDist','dilation'})));

ip.addParamValue('resizeFactor',1,@(x)(numel(x) == 1 & all(x > 0 )));  
ip.addParamValue('logSigma',[2.5],@(x)(numel(x) >0 & all(x > 0 ))); % specify range of nuclei diameters in pixels ie [3 30].
ip.addParamValue('chanRange',[0],@(x)(numel(x) >0 & all(x > 0 ))); %channels for measuring features. If 0, assume all channels.
ip.addParamValue('upSample',2,@(x)(numel(x) == 1 & all(x > 0 )));  
ip.addParamValue('Docker','false',@(x)(ismember(x,{'true','false'})));
ip.addParamValue('dockerParams',0,@(x)(numel(x)==1));

Segmentation label masks for nuclei, cytoplasm, and cell will be saved to a subfolder under each parent image folder as a .tif file. Also saved are a 2-channel tif file with the DAPI and nuclei outlines for quality control. 
