input = "C:/Users/kaabi/Documents/Fiji/Actin/";
output = "C:/Users/kaabi/Documents/Fiji/Output/";


processFolder(input);

setBatchMode(true); 
function processFolder(input) {
	inputDir = File.getName(input);	
	list = getFileList(input);
	for (i = 0; i < list.length; i++) {
			action(input, output, list[i]);
	}
}
// Get Intensity of all channels
function action(input, output, file) {
        open(input + file);
		run("Z Project...", "projection=[Max Intensity]");
		run("Enhance Contrast", "saturated=0.35");
		run("Split Channels");
		
		selectWindow("C4-MAX_"+input+file);
		run("Gaussian Blur...", "sigma=2.5");
		selectWindow("C4-MAX_"+input+file);
		setAutoThreshold("Default dark no-reset");
		run("Analyze Particles...", "size=10-Infinity display exclude clear include summarize overlay add composite");
		
		selectWindow("C3-MAX_"+input+file);
		roiManager("Show All without labels");
		roiManager("Show All with labels");
		roiManager("Measure");
		
		selectWindow("C2-MAX_"+input+file);
		roiManager("Show All without labels");
		roiManager("Show All with labels");
		roiManager("Measure");
		saveAs("Results", output + list[i]+".csv");
		close("*");
} 
//setBatchMode(false);