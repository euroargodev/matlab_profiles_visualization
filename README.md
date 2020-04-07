# MPV: Matlab Profiles Visualization

This application that has been developed by OGS to support the activities related to the **Delayed-Mode Quality Control** (DMQC) of the MOCCA Argo fleet in the Mediterranean and Black Seas.

The application is called **Matlab Profiles Visualization (MPV)** and it is a user friendly tool that reads the Argo NetCDF files, converts the files in MATLAB format, allows Argo profiles selection, produces graphs of temperature and salinity, performs a tailored comparison between the float and reference profiles and provides the main information of float profiles.

**Requirements**

The MPV tool has been developed with MATLAB (version R2017b) and it has been tested on both LINUX and Windows platforms. It requires the use of the MATLAB Mapping Toolbox (https://it.mathworks.com/products/mapping.html). The MPV has been built in MATLAB environment using an Object-Oriented Programming (OOP). The MPV allows to select specific float profiles to be compared with the reference dataset. Moreover, it provides with several diagrams in support of the DMQC analysis.

## TECHNICAL DESCRIPTION

The MPV exploits the OOP features of MATLAB to couple analysis algorithms with a Graphical User Interface (GUI). 

The program differs from common MATLAB scripts since it creates a user-driven space where a user can interact with the app, for instance, by entering text or by clicking a button. Every action is followed by a function call (callback function) and the execution of some code. Callbacks are short functions that must obtain some data to do their job, update the app if necessary and store results for other callbacks. 
The app is just a container and a manager of small functions needed to get a larger task. MATLAB offers a GUIDE, an interface to create interactively the GUI but it provides only a small set of widgets. MPV has been built from scratch by hand and depends on two main classes:
-	**ArgoUI**: it is the main GUI object that creates the main graphical interface and provides the user with every options and information needed to load, visualize and analyse data.
-	**ArgoPlot**: it manages plots of profiles and the interaction with them through the user mouse pointer. For instance, plotted profiles can be selected with the mouse and compared with reference profiles.

The plotting windows are standard MATLAB figure plots. Custom plot windows could be created but they provide less features of the standard MATLAB ones. ArgoUI and ArgoPlot use some external callback functions (defined in a separate source file) to operate on profiles.

## HOW TO USE MPV:

In the "Argo" launcher script, set the following paths:
1. ``root`` (where your reference CTD files are)
1. ``float_data`` (where your float files are)
1. ``dataPath`` (CTD reference files folder)

In your CTD reference files folder, CTD files (in the form: ctd_wmobox.mat) have to be grouped per sub-basins (sub-folders named "Adriatic, Ionian, ...)

In the ``compare`` script the polygon of the main area of interest is loaded ('medPol.mat'):
```matlab
med = load(fullfile(obj.mainApp.paths.data, 'medPol.mat'));
```

The MPV is launched by typing ``argo`` in the MATLAB Command Window

A complete description on how to use the software and examples are given in the report: "D4.4.3 Argo App for data reading, plotting and comparison  v0.2"
