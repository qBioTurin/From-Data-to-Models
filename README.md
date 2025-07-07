# From Data to Models: Analyzing and Integrating Biological Data into Mechanistic Models

Welcome to the GitHub repository for the **PhD Course: From Data to Models**. This course introduces participants to mechanistic modeling in biology, emphasizing how to analyze and integrate real experimental data into computational models. Through lectures and hands-on sessions, students will work with tools such as Petri Nets, epimod, and Flux Balance Analysis (FBA).

[Link Materials](https://drive.google.com/drive/folders/1mGFzbS1nfIQ80TwQxC0FIkY0Mz6y_j7i?usp=drive_link)

---

## Course Structure

The repository is organized into four folders:

- `Day1`: Introduction to mechanistic modeling, Petri Nets, and the epimod framework
- `Day2`: Working with experimental data and defining models
- `Day3`: Model calibration, prediction, and FBA
- `Day4`: Data integration and advanced modeling strategies

Each folder includes:
- Slides
- Scripts
- Hands-on examples and datasets

---

## Course Outline

### Day 1 – Introduction to Mechanistic Modeling and the Importance of Data
- Welcome and course overview
- Introduction to computational models in biology
- Mechanistic modeling foundations: Petri Nets and the epimod framework
- Hands-on session: the Schlögl model, exploring dynamics and sensitivity

### Day 2 – From Case Study to Data Analysis
- From experimental data to model parameters
- Real-world case study: biological context and modeling needs
- Hands-on session: data analysis using ORCA
- Model construction from raw data

### Day 3 – Calibration and Predictive Modeling
- Overview of model calibration techniques
- Performing predictive simulations and what-if analyses
- Introduction to Flux Balance Analysis (FBA)

### Day 4 – Data Integration and Advanced Modeling
- Integrating FBA with experimental data for improved predictions
- Introduction to UnifiedGreatMOD for linking heterogeneous data and models
- Hands-on session: FBA case study with integrated data
- Course wrap-up and summary of key concepts

---

## Required Tools

To run the examples and participate in the hands-on activities, the following tools are required.

Note: Docker and R are already installed on the HPC servers provided for the course. Local installation is optional.

### 1. GreatSPN
Graphical tool for drawing and simulating Petri Nets.  
Installation instructions and downloads are available at:  
[greatspnHOWTOINSTALL](https://github.com/greatspn/SOURCES/blob/master/docs/INSTALL.md)

### 2. Docker
Docker is required to run the `epimod` framework via containerized images.  
Follow the installation guide here: [https://docs.docker.com/get-docker/](https://docs.docker.com/get-docker/)

After installation, enable non-root access to Docker with:

```bash
sudo groupadd docker      # Create the docker group (if not already present)
sudo usermod -aG docker $USER
```

### 3. R, RStudio, and epimod (optional for local execution)
These are not required if you are using the HPC servers provided during the course.

To set up locally:

Install R and RStudio.

Install required packages in R:
```r
install.packages("remotes")
install.packages("fdatest")
remotes::install_github('qBioTurin/epimod', ref='epimod_pFBA', dependencies=TRUE)

library(epimod)
downloadContainers()  # Downloads all necessary Docker images
```

### 4. ORCA

https://github.com/qBioTurin/ORCA
