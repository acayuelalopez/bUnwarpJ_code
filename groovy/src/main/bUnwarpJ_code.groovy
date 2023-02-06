import bunwarpj.BSplineModel
import bunwarpj.MainDialog
import bunwarpj.Mask
import bunwarpj.MiscTools
import bunwarpj.PointHandler
import bunwarpj.Transformation
import ij.IJ
import ij.ImageJ
import ij.ImagePlus
import ij.ImageStack
import ij.measure.ResultsTable
import ij.plugin.ChannelSplitter
import ij.plugin.RGBStackMerge
import ij.process.ImageProcessor
import ij.process.LUT
import loci.plugins.BF
import loci.plugins.in.ImporterOptions
import java.awt.Point


// INPUT UI
//
#@File(label="Input File Directory", style="directory") inputFile
#@File(label="Beads reference File", style="file") beadsFile
#@File(label = "Output directory", style="directory") outputDir
#@Integer(label="Fixed Channel", value=1) fixedCh
#@Boolean (label="Run headless", default="false") headless
//Advanced Parameters
#@Integer(label="Registration Mode", value=1) mode
#@Integer(label="Image Subsample Factor", value=0) img_subsamp_fact
#@Integer(label="Initial Deformation", value=0) min_scale_deformation
#@Integer(label="Final Deformation", value=2) max_scale_deformation
#@Double(label="Divergence Weight", value=0.0) divWeight
#@Double(label="Curl Weight", value=0.0) curlWeight
#@Double(label="Landmark Weight", value=0.0) landmarkWeight
#@Double(label="Image Weight", value=1.0) imageWeight
#@Double(label="Consistency Weight", value=10.0) consistencyWeight
#@Double(label="Stop Threshold", value=0.01) stopThreshold

// IDE
//
//def inputFile = new File("/home/anaacayuela/Ana_pruebas_imageJ/bUnwarpJ/paper/set_2/images/Sylvia Gutierrez-Erlandsson - Hoechst 1x_34ºC.lif");
//def beadsFile = new File("/home/anaacayuela/Ana_pruebas_imageJ/bUnwarpJ/paper/set_2/ref_beads/Sylvia Gutierrez-Erlandsson - Beads ref.lif");
//def outputDir = new File("/home/anaacayuela/Ana_pruebas_imageJ/bUnwarpJ/paper/set_2/results");
//def fixedCh = 0;
//def headless = true;
//new ImageJ().setVisible(false);


//bUnwarpJ Advanced parameters
//int mode = 1;
//int img_subsamp_fact = 0;
//int min_scale_deformation = 0;
//int max_scale_deformation = 2
//def divWeight = 0.0.doubleValue();
//def curlWeight = 0.0.doubleValue();
//def landmarkWeight = 0.0.doubleValue();
//def imageWeight = 1.0.doubleValue();
//def consistencyWeight = 10.0.doubleValue();
//def stopThreshold = 0.01.doubleValue();

//Get image file to be analyzed

if (inputFile.getName().endsWith(".tif")) {
    //Declare each image to process within input directory
    def imp = new ImagePlus(inputFile.getAbsolutePath());
    def impTitle = null;
    if (imp.getTitle().contains("/")) {
        impTitle = imp.getTitle().replaceAll("/", "").replaceAll(".tif", "");
    } else {
        impTitle = imp.getTitle().replaceAll(".tif", "")
    }
    //Define output directory per image
    def outputImageDir = new File(
            outputDir.getAbsolutePath() + File.separator + impTitle);

    if (!outputImageDir.exists()) {
        def results = false;

        try {
            outputImageDir.mkdir();
            results = true;
        } catch (SecurityException se) {
            // handle it
        }
    }
    //Declare arrays to save moving images, channel and min/max display range
    def movingImp = new ArrayList<ImagePlus>();
    def movingImpCh = new ArrayList<String>();
    def minMovingImp = new ArrayList<Integer>();
    def maxMovingImp = new ArrayList<Integer>();
    def lutMovingImp = new ArrayList<LUT>();

    //Get calibration from non-transformed image
    def cal = imp.getCalibration();
    //Split channels
    def channels = ChannelSplitter.split(imp);
    //Create array to save each channel (including transformed)
    def channelsToMerge = new ImagePlus[channels.length];
    //Declare fixed image (target image)
    def fixedImp = channels[fixedCh].duplicate();
    IJ.resetMinAndMax(fixedImp)
    //def lutFixedImp = fixedImp.getProcessor().getLut();
    //Save fixed image channel
    IJ.saveAs(fixedImp, "Tiff", outputImageDir.getAbsolutePath()
            + File.separator + "Fixed_Ch_" + fixedCh.toString() + "_" + impTitle);
    //Add fixed channel to array of images to merge
    channelsToMerge[fixedCh] = fixedImp
    //Iterate through channels
    for (def channel = 0.intValue(); channel < channels.length; channel++) {
        //Iterate to save moving channels except from fixed channels
        if (channel !== fixedCh) {
            IJ.resetMinAndMax(channels[channel]);
            movingImp.add(channels[channel].duplicate());
            movingImpCh.add(channel.toString());
            minMovingImp.add(channels[channel].getDisplayRangeMin().toInteger());
            maxMovingImp.add(channels[channel].getDisplayRangeMax().toInteger());
            lutMovingImp.add(channels[channel].getProcessor().getLut())
        }
    }

    //Create beads reference image
    def options = new ImporterOptions();
    options.setId(beadsFile.getAbsolutePath());
    options.setSplitChannels(false);
    options.setSplitTimepoints(false);
    options.setSplitFocalPlanes(false);
    options.setAutoscale(false);
    options.setStackFormat(ImporterOptions.VIEW_HYPERSTACK);
    options.setStackOrder(ImporterOptions.ORDER_XYCZT);
    options.setCrop(false);
    options.setOpenAllSeries(true);
    def impsBeads = BF.openImagePlus(options);
    //Get first serie of beads file as reference beads file
    def impBeads = impsBeads[0];
    //Split channels
    def beadsChannels = ChannelSplitter.split(impBeads);
    //Define target image from beads
    def targetBeads = beadsChannels[fixedCh.intValue()].duplicate();
    IJ.resetMinAndMax(targetBeads)
    //Define source image from beads
    def sourceCh = 0.intValue();
    if (fixedCh == 0) {
        sourceCh = 1.intValue();
    } else {
        sourceCh = 0.intValue();
    }

    def sourceBeads = beadsChannels[sourceCh].duplicate();
    IJ.resetMinAndMax(sourceBeads)
    //Get transformation object that contains all the registration information.
    def transf = computeTransformationBatch(targetBeads, sourceBeads, targetBeads.getMask(), sourceBeads.getMask(), mode, img_subsamp_fact, min_scale_deformation, max_scale_deformation, divWeight, curlWeight, landmarkWeight, imageWeight, consistencyWeight, stopThreshold);
    //Get number of intervals from tranformation object between B-spline coefficients
    def intervals = transf.getIntervals();
    //Get the direct deformation X coefficients
    def cx = transf.getDirectDeformationCoefficientsX()
    //Get the direct deformation Y coefficients.
    def cy = transf.getDirectDeformationCoefficientsY()
    //Iterate through moving channels
    for (int j = 0; j < movingImp.size(); j++) {
        //Create output stack per moving channel
        def outputStackMoving = new ImageStack(movingImp.get(j).getWidth(),
                movingImp.get(j).getHeight());
        // apply transform to each slice of the stack
        for (int z = 1; z <= movingImp.get(j).getImageStackSize(); z++) {
            def ip = movingImp.get(j).getImageStack().getProcessor(z);

            def source = new BSplineModel(ip, false, 1);

            def movingImage = new ImagePlus("", ip);

            def result = MiscTools.applyTransformationMT(movingImage, fixedImp,
                    source, intervals, cx, cy);

            outputStackMoving.addSlice("", result);
        }
        //Save each moving stack as ImagePlus
        def impFinalMoving = new ImagePlus("", outputStackMoving);
        impFinalMoving.setDisplayRange(minMovingImp.get(j).intValue(), maxMovingImp.get(j).intValue());
        impFinalMoving.setCalibration(cal);
        impFinalMoving.setLut(lutMovingImp.get(j));
        IJ.run(impFinalMoving, String.format("%s-bit", fixedImp.getBitDepth().toString()), "");
        //Save moving image channel
        IJ.saveAs(impFinalMoving, "Tiff", outputImageDir.getAbsolutePath()
                + File.separator + "Moving_Ch_" + movingImpCh.get(j) + "_" + impTitle);
        //Add moving channel to array to merge channels
        channelsToMerge[movingImpCh.get(j).toInteger().intValue()] = impFinalMoving;
    }
    //Merge moving and fixed channels as ImagePlus
    def impMerged = new RGBStackMerge().mergeHyperstacks(channelsToMerge, false);
    impMerged.setCalibration(cal);
    IJ.resetMinAndMax(impMerged);
    //Save merged channels image (aligned)
    IJ.saveAs(impMerged, "Tiff", outputImageDir.getAbsolutePath()
            + File.separator + "Aligned_" + impTitle);

}
if (inputFile.getName().endsWith(".lif")) {
    def options = new ImporterOptions();
    options.setId(inputFile.getAbsolutePath());
    options.setSplitChannels(false);
    options.setSplitTimepoints(false);
    options.setSplitFocalPlanes(false);
    options.setAutoscale(false);
    options.setStackFormat(ImporterOptions.VIEW_HYPERSTACK);
    options.setStackOrder(ImporterOptions.ORDER_XYCZT);
    options.setCrop(false);
    options.setOpenAllSeries(true);
    def imps = BF.openImagePlus(options);

    for (int j = 0; j < imps.length; j++) {
        //Declare each image to process within input directory
        def imp = imps[j];
        def impTitleSerie = null;
        if (imp.getTitle().contains("/")) {
            impTitleSerie = imp.getTitle().replaceAll("/", "").replaceAll(".lif", "") + "_serie_" + (j + 1).toString();
        } else {
            impTitleSerie = imp.getTitle().replaceAll(".lif", "") + "_serie_" + (j + 1).toString();
        }

        //Define output directory per image
        def outputImageParentDir = new File(
                outputDir.getAbsolutePath() + File.separator + inputFile.getName().replaceAll(".lif", ""));

        if (!outputImageParentDir.exists()) {
            def results = false;

            try {
                outputImageParentDir.mkdir();
                results = true;
            } catch (SecurityException se) {
                // handle it
            }
        }
        def outputImageDir = new File(
                outputImageParentDir.getAbsolutePath() + File.separator + impTitleSerie);

        if (!outputImageDir.exists()) {
            def results = false;

            try {
                outputImageDir.mkdir();
                results = true;
            } catch (SecurityException se) {
                // handle it
            }
        }
        //Declare arrays to save moving images, channel and min/max display range
        def movingImp = new ArrayList<ImagePlus>();
        def movingImpCh = new ArrayList<String>();
        def minMovingImp = new ArrayList<Integer>();
        def maxMovingImp = new ArrayList<Integer>();
        def lutMovingImp = new ArrayList<LUT>();

        //Get calibration from non-transformed image
        def cal = imp.getCalibration();
        //Split channels
        def channels = ChannelSplitter.split(imp);
        //Create array to save each channel (including transformed)
        def channelsToMerge = new ImagePlus[channels.length];
        //Declare fixed image (target image)
        def fixedImp = channels[fixedCh].duplicate();
        IJ.resetMinAndMax(fixedImp)
        //Save fixed image channel
        IJ.saveAs(fixedImp, "Tiff", outputImageDir.getAbsolutePath()
                + File.separator + "Fixed_Ch_" + fixedCh.toString() + "_" + impTitleSerie);
        //Add fixed channel to array of images to merge
        channelsToMerge[fixedCh] = fixedImp
        //Iterate through channels
        for (def channel = 0.intValue(); channel < channels.length; channel++) {
            //Iterate to save moving channels except from fixed channels
            if (channel !== fixedCh) {
                IJ.resetMinAndMax(channels[channel]);
                movingImp.add(channels[channel].duplicate());
                movingImpCh.add(channel.toString());
                minMovingImp.add(channels[channel].getDisplayRangeMin().toInteger());
                maxMovingImp.add(channels[channel].getDisplayRangeMax().toInteger());
                lutMovingImp.add(channels[channel].getProcessor().getLut())
            }
        }
        //Create beads reference image
        options = new ImporterOptions();
        options.setId(beadsFile.getAbsolutePath());
        options.setSplitChannels(false);
        options.setSplitTimepoints(false);
        options.setSplitFocalPlanes(false);
        options.setAutoscale(false);
        options.setStackFormat(ImporterOptions.VIEW_HYPERSTACK);
        options.setStackOrder(ImporterOptions.ORDER_XYCZT);
        options.setCrop(false);
        options.setOpenAllSeries(true);
        def impsBeads = BF.openImagePlus(options);
        //Get first serie of beads file as reference beads file
        def impBeads = impsBeads[0];
        //Split channels
        def beadsChannels = ChannelSplitter.split(impBeads);
        //Define target image from beads
        def targetBeads = beadsChannels[fixedCh.intValue()].duplicate();
        IJ.resetMinAndMax(targetBeads)
        //Define source image from beads
        def sourceCh = 0.intValue();
        if (fixedCh == 0) {
            sourceCh = 1.intValue();
        } else {
            sourceCh = 0.intValue();
        }

        def sourceBeads = beadsChannels[sourceCh].duplicate();
        IJ.resetMinAndMax(sourceBeads)
        //Get transformation object that contains all the registration information.
        def transf = computeTransformationBatch(targetBeads, sourceBeads, targetBeads.getMask(), sourceBeads.getMask(), mode, img_subsamp_fact, min_scale_deformation, max_scale_deformation, divWeight, curlWeight, landmarkWeight, imageWeight, consistencyWeight, stopThreshold);
        //Get number of intervals from tranformation object between B-spline coefficients
        def intervals = transf.getIntervals();
        //Get the direct deformation X coefficients
        def cx = transf.getDirectDeformationCoefficientsX()
        //Get the direct deformation Y coefficients.
        def cy = transf.getDirectDeformationCoefficientsY()
        //Iterate through moving channels
        for (int z = 0; z < movingImp.size(); z++) {
            //Create output stack per moving channel
            def outputStackMoving = new ImageStack(movingImp.get(z).getWidth(),
                    movingImp.get(z).getHeight());
            // apply transform to each slice of the stack
            for (int x = 1; x <= movingImp.get(z).getImageStackSize(); x++) {
                def ip = movingImp.get(z).getImageStack().getProcessor(x);

                def source = new BSplineModel(ip, false, 1);

                def movingImage = new ImagePlus("", ip);

                def result = MiscTools.applyTransformationMT(movingImage, fixedImp,
                        source, intervals, cx, cy);

                outputStackMoving.addSlice("", result);
            }
            //Save each moving stack as ImagePlus
            def impFinalMoving = new ImagePlus("", outputStackMoving);
            impFinalMoving.setDisplayRange(minMovingImp.get(z).intValue(), maxMovingImp.get(z).intValue());
            impFinalMoving.setCalibration(cal);
            IJ.run(impFinalMoving, String.format("%s-bit", fixedImp.getBitDepth().toString()), "");
            impFinalMoving.setLut(lutMovingImp.get(z));
            //Save moving image channel
            IJ.saveAs(impFinalMoving, "Tiff", outputImageDir.getAbsolutePath()
                    + File.separator + "Moving_Ch_" + movingImpCh.get(z) + "_" + impTitleSerie);
            //Add moving channel to array to merge channels
            channelsToMerge[movingImpCh.get(z).toInteger().intValue()] = impFinalMoving;
        }
        //Merge moving and fixed channels as ImagePlus
        def impMerged = RGBStackMerge.mergeChannels(channelsToMerge, false);
        impMerged.setCalibration(cal);
        IJ.resetMinAndMax(impMerged);
        //Save merged channels image (aligned)
        IJ.saveAs(impMerged, "Tiff", outputImageDir.getAbsolutePath()
                + File.separator + "Aligned_" + impTitleSerie);

    }
}


// exit
//
if (headless)
    System.exit(0)

Transformation computeTransformationBatch(ImagePlus targetImp,
                                          ImagePlus sourceImp,
                                          ImageProcessor targetMskIP,
                                          ImageProcessor sourceMskIP,
                                          int mode,
                                          int img_subsamp_fact,
                                          int min_scale_deformation,
                                          int max_scale_deformation,
                                          double divWeight,
                                          double curlWeight,
                                          double landmarkWeight,
                                          double imageWeight,
                                          double consistencyWeight,
                                          double stopThreshold) {
    // Produce side information
    def imagePyramidDepth = max_scale_deformation - min_scale_deformation + 1;
    def min_scale_image = 0.intValue();

    // output level to -1 so nothing is displayed
    def outputLevel = -1.intValue();

    def showMarquardtOptim = false;

    // Create target image model
    def target = new BSplineModel(targetImp.getProcessor(), true,
            Math.pow(2, img_subsamp_fact).intValue());

    target.setPyramidDepth(imagePyramidDepth + min_scale_image);
    target.startPyramids();

    // Create target mask
    def targetMsk = (targetMskIP != null) ? new Mask(targetMskIP, true)
            : new Mask(targetImp.getProcessor(), false);

    def targetPh = null;

    // Create source image model
    def bIsReverse = true;

    def source = new BSplineModel(sourceImp.getProcessor(), bIsReverse,
            Math.pow(2, img_subsamp_fact).intValue());

    source.setPyramidDepth(imagePyramidDepth + min_scale_image);
    source.startPyramids();

    // Create source mask
    def sourceMsk = (sourceMskIP != null) ? new Mask(sourceMskIP, true)
            : new Mask(sourceImp.getProcessor(), false);

    def sourcePh = null;

    // Load points rois as landmarks if any.
    def sourceStack = new Stack<Point>();
    def targetStack = new Stack<Point>();
    MiscTools.loadPointRoiAsLandmarks(sourceImp, targetImp, sourceStack, targetStack);

    sourcePh = new PointHandler(sourceImp);
    targetPh = new PointHandler(targetImp);

    while ((!sourceStack.empty()) && (!targetStack.empty())) {
        def sourcePoint = (Point) sourceStack.pop();
        def targetPoint = (Point) targetStack.pop();
        sourcePh.addPoint(sourcePoint.x, sourcePoint.y);
        targetPh.addPoint(targetPoint.x, targetPoint.y);
    }


    // Set no initial affine matrices
    def sourceAffineMatrix = null;
    def targetAffineMatrix = null;

    // Join threads
    try {
        source.getThread().join();
        target.getThread().join();
    }
    catch (InterruptedException e) {
        IJ.error("Unexpected interruption exception " + e);
    }

    // Perform registration
    def output_ip = new ImagePlus[2];
    output_ip[0] = null;
    output_ip[1] = null;

    MainDialog dialog = null;
    def originalSourceIP = sourceImp.getProcessor();
    def originalTargetIP = targetImp.getProcessor();

    // Setup registration parameters
    def warp = new Transformation(
            sourceImp, targetImp, source, target, sourcePh, targetPh,
            sourceMsk, targetMsk, sourceAffineMatrix, targetAffineMatrix,
            min_scale_deformation, max_scale_deformation, min_scale_image, divWeight,
            curlWeight, landmarkWeight, imageWeight, consistencyWeight, stopThreshold,
            outputLevel, showMarquardtOptim, mode, null, null, output_ip[0], output_ip[1], dialog,
            originalSourceIP, originalTargetIP);
    if (mode == 1.intValue())
        warp.doUnidirectionalRegistration();
    else
        warp.doBidirectionalRegistration();


    return warp;

} // end computeTransformationBatch
