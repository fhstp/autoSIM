## autoSIM: a Matlab-based fully automated workflow for semi-personalized physics-based OpenSim simulations.

This is a bundle of Matlab scripts and files that allow to run different *.osim models on a set of *.c3d files. It will help to handle all of the necessary files (*.xml, *.trc, *.mot, etc.) to run the simulations. Once set up, it is possible to iterate several trials in Matlab automatically. AutoSim was originally developed for the Joint Articulation Mechanics (JAM, https://github.com/clnsmith/opensim-jam) workflow of C. Smith, but then was extended to also run other models.

AutoSim was developed to handle large-scale datasets and to fully automate OpenSim simulations across vast amounts of data, addressing the growing demand for big data in the era of machine learning and biomechanics. For example, we successfully used AutoSim to generate semi-personalized simulations of ~80,000 gait cycles using the "Lernergopal model", all within less than two weeks of processing time.

## Prerequisites

Before you begin, ensure you have met the following requirements:

- You will need to have OpenSim 4.0 installed and the Matlab scripting environment using the API set up, see https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scripting+with+Matlab 
- Code was tested with Matlab R2024a and WIN11 on Intel machines.

## In Matlab you will need the following toolboxes:
- Image Processing Toolbox — Version 11.6+ (R2022b)  
- Parallel Computing Toolbox — Version 7.7+ (R2022b)  
- Phased Array System Toolbox — Version 4.8+ (R2022b)  
- Signal Processing Toolbox — Version 9.1+ (R2022b)  
- Statistics and Machine Learning Toolbox — Version 12.4+ (R2022b)  

## Clone the Repository

To get started, clone the repository with:

```bash
git clone https://github.com/bhorsak/autoSIM.git
```

## Using autoSIM

- The code was primarily developed for gait data using events. For now autoSIM expects events to have these lables: 'R_Foot_Strike' or 'L_Foot_Strike' and 'R_Foot_Off' or 'L_Foot_Off' (see `.\_commonFiles\sharedFunctions\getEventsOSS.m`)
- To test the workflow you can use the example data provided in `.\_commonFiles\dataExamples`.
- Each model workflow has its own `start.m` file. In each file, there are some hardcoded settings which you will have to set the first time before running autoSIM, e.g. marker set definition. 
- The entire repo expects a Vicon Nexus database folder. To run autoSIM your folder at least needs to contain *.c3d files for the condition to analyze and a static file. All *.c3d files need to include the condition information within each filename, such as `walking` (e.g., 001_walking_file.c3d). You will also need the `enf.files` which hold the information of valid/invalid force plate hits and which foot hit which force plate. If you don't have such a file you may want to create it yourself. The `start.m` files assumes that all necessary *.c3d files are located in one single folder (e.g. of one subject). It finds all relevant *.c3d files, creates the input data (e.g. start and stop time, body side, etc.), creates the necessary *.trc and *.mot files, and the external loads file and runs the workflow for each file. Note that, when a *.c3d file holds two valid force plate contacts (e.g. left and right). Simulations are performed separately for each contact.
- In case you want to make use of the built-in personalization methods, such as the Torsion-Tool a `data.xml` for each session folder is needed. It contains in a standardized format information such as the tibio-femoral alignment angle in the frontal plane or tibial torsion. For detaisl please review the example file in `.\_commonFiles\dataExamples`.
- The repository also includes additional scripts such as `postprocessResults.m` which runs a set of post processing steps. For details see the header of each file. 

## Adding New Markerset and Gait Lab Configurations

Currently, this process is not fully automated. Users must manually add their configuration in specific code sections and corresponding folders containing setup scripts and hardcoded settings.

**Part 1: Setup Files**  
Locate the `setupFiles` folder within either the `comak` or `varModels` branch. Inside, you will find subfolders with template setup files for different gait labs and markersets, for example: `\comak\setupFiles\Models\OSS_FHSTP`


Use these as templates for your new markerset. Explore the contained setup files, which should be familiar to OpenSim users. The easiest approach is to copy an existing folder that works with the example data (e.g., in `.\_commonFiles\dataExamples\OSS`), then rename it to match your lab.
To avoid confusion, rename key files like the markerset file, for instance to `MyLab_Cleveland_MarkerSet_COMAK.xml`


Adapt the setup files accordingly, paying particular attention to:  
- The markerset file  
- The scaling template file  
- The inverse kinematics settings file (make sure your new markerset is hardcoded into the setup file)  

Once adapted, Part 1 is complete.

**Part 2: Force Plate and Marker Trajectory Configuration**  
The function `._commonFiles\sharedFunctions\c3d2OpenSim.m` handles creation of external loads and `.mot` files containing force plate and `.trc` marker data. You need to define your lab within this function. AutoSIM currently uses `switch` statements to handle different labs. Add your lab with a unique acronym (e.g., `"MyLab"`) and adjust the rotation matrices to align your data correctly. A practical method: set a breakpoint right after calling `c3d2OpenSim`, define your lab, make adjustmenets to the rotation matrices, run the function, then verify in OpenSim by loading a *.mot and *.trc file of a trial to check if ground reaction forces and marker data are correctly aligned.

**Part 3: Integration of your lab into autoSIM**  
This step can be tricky but mainly involves adding switch statements with your lab-tag (e.g. "MyLab") and filenames. Enable "Pause on Error" in MATLAB to identify functions where your lab is not yet fully integrated (mainly related to missing `switch` cases).

Files you may need to adjust include (in order of being called in autoSIM):  
- `startComak.m` or `startVarModels.m`: Look for the "HARDCODED Settings" section near the end.  
- `prepareInputData.m`: Adjust the `switch` statement around line 269 for "Get events based on lab data".  
- `prepareScaledModel.m`: Modify the top-level `switch` statement.  
- `run_secConstrainSimulation` (comak only): Update the top `switch` statement.  
- `osimjam_workflow.m` or `Model_workflow.m`: Modify the top-level `switch` statement.

## Contributors

Thanks to the following people who gave me incredible supported for developing this repository:

- [@BryceKillian](https://bitbucket.org/BKillen/)
- [@IlseJonkers](https://www.kuleuven.be/wieiswie/en/person/00015567)
- [@SamVanRossom](https://www.kuleuven.be/wieiswie/en/person/00093377)
- [@MarkSimonlehner](https://www.fhstp.ac.at/de/uber-uns/mitarbeiter-innen-a-z/simonlehner-mark)
- [@BernhardDumphart](https://www.fhstp.ac.at/de/uber-uns/mitarbeiter-innen-a-z/dumphart-bernhard)
- [@BernhardGuggenberger](bernhard.guggenberger2@fh-joanneum.at)
- [@HansKainz](hans.kainz@univie.ac.at)
- [@WilliKoller](willi.koller@univie.ac.az)

## Citation
Please cite the following reference(s) if you use autoSIM:

Horsak, B., Krondorfer, P., Unglaube, F., Slijepčević, D., Kranzl, A., 2025. Feasibility of fully automated and semi-personalized musculoskeletal simulations to process large-scale gait datasets, in: Gait & Posture, GAMMA Congress 2025 (26 – 29 March 2025, St Gallen Switzerland). pp. S18–S19. https://doi.org/10.1016/j.gaitpost.2025.01.059

```
@article{HORSAK2025S18,
title = {Feasibility of fully automated and semi-personalized musculoskeletal simulations to process large-scale gait datasets},
journal = {Gait & Posture},
volume = {117},
pages = {S18-S19},
year = {2025},
note = {GAMMA Congress 2025 (26 – 29 March 2025, St Gallen Switzerland)},
issn = {0966-6362},
doi = {https://doi.org/10.1016/j.gaitpost.2025.01.059},
url = {https://www.sciencedirect.com/science/article/pii/S0966636225000591},
author = {Brian Horsak and Philipp Krondorfer and Fabian Unglaube and Djordje Slijepčević and Andreas Kranzl}
}
```

Koller, W., Horsak, B., Kranzl, A., Unglaube, F., Baca, A., Kainz, H., 2025. Physiological plausible muscle paths: A MATLAB tool for detecting and resolving muscle path discontinuities in musculoskeletal OpenSim models, in: Gait & Posture, GAMMA Congress 2025 (26 – 29 March 2025, St Gallen Switzerland). pp. S21–S22. https://doi.org/10.1016/j.gaitpost.2025.01.063

```
@inproceedings{koller_2025_physiological,
  title = {Physiological Plausible Muscle Paths: {{A MATLAB}} Tool for Detecting and Resolving Muscle Path Discontinuities in Musculoskeletal {{OpenSim}} Models},
  shorttitle = {Physiological Plausible Muscle Paths},
  booktitle = {Gait \& {{Posture}}},
  author = {Koller, Willi and Horsak, Brian and Kranzl, Andreas and Unglaube, Fabian and Baca, Arnold and Kainz, Hans},
  year = {2025},
  series = {{{GAMMA Congress}} 2025 (26 -- 29 {{March}} 2025, {{St Gallen Switzerland}})},
  volume = {117},
  pages = {S21-S22},
  doi = {10.1016/j.gaitpost.2025.01.063}
}
```

## Contact
If you want to contact me you can reach me at <brian.horsak@fhstp.ac.at>.

## License
This project is released under the [MIT License]. Copyright (c) 2024 Brian Horsak Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. 

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 

## Third-Party Code and Licenses,
see (`.\_commonFiles\sharedFunctions\ThirdPartyFiles`) and included OpenSim models. 

This project also includes or depends on third-party components with their own licenses (e.g., BSD, MPL, CC-BY-NC). These components are clearly marked in the repository. Users of this project are responsible for ensuring compliance with the corresponding licenses when using or redistributing those parts.

## Further References
Depending on which Model or built-in Tool you will be using, please also consider to cite the following papers. Note the list might not be complete. Please double-check the settings you make and make sure to cite Third-Party Tools, Code, Models, ... appropriately.

```
@article{laiWhyAreAntagonist2017,
  title = {Why Are Antagonist Muscles Co-Activated in My Simulation? {{A}} Musculoskeletal Model for Analysing Human Locomotor Tasks},
  shorttitle = {Why Are Antagonist Muscles Co-Activated in My Simulation?},
  author = {Lai, Adrian K.M. and Arnold, Allison S. and Wakeling, James M.},
  year = {2017},
  journal = {Annals of biomedical engineering},
  volume = {45},
  number = {12},
  pages = {2762--2774},
  issn = {0090-6964},
  doi = {10.1007/s10439-017-1920-7},
  pmcid = {PMC5989715},
  pmid = {28900782}
}

@article{lenhartPredictionValidationLoadDependent2015,
  title = {Prediction and {{Validation}} of {{Load-Dependent Behavior}} of the {{Tibiofemoral}} and {{Patellofemoral Joints During Movement}}},
  author = {Lenhart, Rachel L. and Kaiser, Jarred and Smith, Colin R. and Thelen, Darryl G.},
  year = {2015},
  journal = {Annals of Biomedical Engineering},
  volume = {43},
  number = {11},
  pages = {2675--2685},
  issn = {1573-9686},
  doi = {10.1007/s10439-015-1326-3},
  langid = {english}
}

@article{lernerHowTibiofemoralAlignment2015,
  title = {How Tibiofemoral Alignment and Contact Locations Affect Predictions of Medial and Lateral Tibiofemoral Contact Forces},
  author = {Lerner, Zachary F. and DeMers, Matthew S. and Delp, Scott L. and Browning, Raymond C.},
  year = {2015},
  journal = {Journal of Biomechanics},
  volume = {48},
  number = {4},
  pages = {644--650},
  issn = {00219290},
  doi = {10.1016/j.jbiomech.2014.12.049},
  langid = {english}
}

@article{rajagopalFullBodyMusculoskeletalModel2016,
  title = {Full-{{Body Musculoskeletal Model}} for {{Muscle-Driven Simulation}} of {{Human Gait}}},
  author = {Rajagopal, Apoorva and Dembia, Christopher L. and DeMers, Matthew S. and Delp, Denny D. and Hicks, Jennifer L. and Delp, Scott L.},
  year = {2016},
  journal = {IEEE transactions on bio-medical engineering},
  volume = {63},
  number = {10},
  pages = {2068--2079},
  issn = {1558-2531},
  doi = {10.1109/TBME.2016.2586891},
  langid = {english},
  pmcid = {PMC5507211},
  pmid = {27392337}
}

@article{smithCanAlteredNeuromuscular2019,
  title = {Can Altered Neuromuscular Coordination Restore Soft Tissue Loading Patterns in Anterior Cruciate Ligament and Menisci Deficient Knees during Walking?},
  author = {Smith, Colin R. and Brandon, Scott C. E. and Thelen, Darryl G.},
  year = {2019},
  journal = {Journal of Biomechanics},
  volume = {82},
  pages = {124--133},
  issn = {0021-9290},
  doi = {10/gmn54k},
  langid = {english}
}

@article{smithInfluenceComponentAlignment2016,
  title = {The {{Influence}} of {{Component Alignment}} and {{Ligament Properties}} on {{Tibiofemoral Contact Forces}} in {{Total Knee Replacement}}},
  author = {Smith, Colin R. and Vignos, Michael F. and Lenhart, Rachel L. and Kaiser, Jarred and Thelen, Darryl G.},
  year = {2016},
  journal = {Journal of Biomechanical Engineering},
  volume = {138},
  number = {2},
  pages = {021017},
  issn = {0148-0731, 1528-8951},
  doi = {10.1115/1.4032464},
  langid = {english}
}

@article{uhlrichMuscleCoordinationRetraining2022a,
  title = {Muscle Coordination Retraining Inspired by Musculoskeletal Simulations Reduces Knee Contact Force},
  author = {Uhlrich, Scott D. and Jackson, Rachel W. and Seth, Ajay and Kolesar, Julie A. and Delp, Scott L.},
  year = {2022},
  journal = {Scientific Reports},
  volume = {12},
  number = {1},
  pages = {9842},
  publisher = {Nature Publishing Group},
  issn = {2045-2322},
  doi = {10.1038/s41598-022-13386-9},
  langid = {english}
}

@article{veerkampTorsionToolAutomated2021,
  title = {Torsion {{Tool}}: {{An}} Automated Tool for Personalising Femoral and Tibial Geometries in {{OpenSim}} Musculoskeletal Models},
  shorttitle = {Torsion {{Tool}}},
  author = {Veerkamp, Kirsten and Kainz, Hans and Killen, Bryce A. and J{\'o}nasd{\'o}ttir, Hulda and {van der Krogt}, Marjolein M.},
  year = {2021},
  journal = {Journal of Biomechanics},
  volume = {125},
  pages = {110589},
  issn = {0021-9290},
  doi = {10.1016/j.jbiomech.2021.110589},
  langid = {english}
}
```





