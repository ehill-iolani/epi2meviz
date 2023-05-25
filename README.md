# EPI2MEVIZ

A R based shiny app that is intendted to be used as a companion analysis to the 16S EPI2ME workflow from Oxford Nanopore Technologies.

This app will take the .csv file downloaded from the EPI2ME analysis and will conduct 3 secondary analyses:
1) Generate rarefaction curves to determine if the run has achieved the proper sampling depth
2) Calculate the relative abundance of each species detected in the EPI2ME results and generate the corresponding plot
3) Conduct a Bray Curtis dissimilarity principle coordinate analysis (PCoA) and plot PCoA 1 vs PCoA 2 (only works if there are > 2 barcodes in the analysis)

You will need to upload the .csv file downloaded from the EPI2ME analysis, set an accuracy threshold based on the EPI2ME results, and select the relevant barcodes of interest.

# Installation

Pull the Docker image from Dockerhub to run this pipeline as a container:
```
docker pull ethill/epi2meviz:alpha
docker run --name=epi2meviz --rm -d -p 3838:3838 ethill/epi2meviz:alpha
```

Open a browser and nativate to http://localhost:3838/ to use the application.

Alterantively, there is a dockerfile included in the repository which allows you to build the image yourself.
Navigate to where you would store the github repository and then enter the following:
```
git clone https://github.com/ehill-iolani/epi2meviz.git
cd epi2meviz
docker build -t epi2meviz:local .
docker run --name=epi2meviz --rm -d -p 3838:3838 epi2meviz:local
```

Open a browser and nativate to http://localhost:3838 to use the application.
