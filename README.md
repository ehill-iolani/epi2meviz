# EPI2MEVIZ

A R based shiny app that is intended to be used as a companion analysis to the 16S EPI2ME workflow from Oxford Nanopore Technologies.

This app will take the .csv file downloaded from the EPI2ME analysis and will conduct 3 secondary analyses:
1) Generate rarefaction curves to determine if the run has achieved the proper sampling depth
2) Calculate the relative abundance of each species detected in the EPI2ME results and generate the corresponding plot
3) Conduct a Bray Curtis dissimilarity principle coordinate analysis (PCoA) and plot PCoA 1 vs PCoA 2 (only works if there are > 2 barcodes in the analysis)

You will need to upload the .csv file downloaded from the EPI2ME analysis, set an accuracy threshold based on the EPI2ME results, and select the relevant barcodes of interest.

# Usage

To use this app you must have access to the .csv file generated from the EPI2ME analysis and access to the EPI2ME dashboard for the analysis stats.
Once you have access to both do the following:
1) Upload the .csv to the app using the "Browse" button
2) Set the "Average EPI2ME Accuracy" value to the average accuracy reported in the EPI2ME dashboard
3) Select the barcodes in you sequencing runs from the checkboxes under "Barcodes to Analyze"
4) Click the "submit" button

Once the analysis completes, the tabs at the top of the page will show the plots/tables generated from the analysis.

You can then download the .pdf of the plot of interest by clicking the "Download PDF" button in the upper right.

# Local Installation

Pull the Docker image from Dockerhub to run this pipeline as a container:
```
docker pull ethill/epi2meviz:beta
docker run --name=epi2meviz --rm -d -p 3838:3838 ethill/epi2meviz:beta
```

Open a browser and nativate to http://localhost:3838/ to use the application.

Alterantively, there is a dockerfile included in the repository which allows you to build the image yourself.
Navigate to where you will store the github repository and then enter the following:
```
git clone https://github.com/ehill-iolani/epi2meviz.git
cd epi2meviz
docker build -t epi2meviz:local .
docker run --name=epi2meviz --rm -d -p 3838:3838 epi2meviz:local
```

Open a browser and nativate to http://localhost:3838/ to use the application.
