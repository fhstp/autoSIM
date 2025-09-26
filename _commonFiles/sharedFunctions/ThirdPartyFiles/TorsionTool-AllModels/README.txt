%                                                                         %
%    Authors: Hulda Jónasdóttir & Kirsten Veerkamp                        %
%                            February 2021                                %
%    email:    k.veerkamp@amsterdamumc.nl                                 % 
%                                                                         %
%  modified by Elias Wallnöfer & Willi Koller   (Mai 2023)                %
%                                                                         %
%  email: willi.koller@univie.ac.at                                       %
%                                                                         %

1. Open the main script (MAIN_TorsionToolAllModels).
2. Give the subject-specific femoral anteversion (AV) and neck-shaft (NS) angles,
	as well as the tibial torsion (TT) angles, as input for the right and left leg.
	Lines which require these inputs are indicated by a % at the end of the line.
3. The final model with personalised torsions is saved in the DEFORMED_MODEL
	folder, and is called FINAL_PERSONALISEDTORSIONS.osim.
	The adjusted markerset can also be found in this folder.

note1: The angle definitions for AV and TT are as follows:
	- AV: positive: femoral anteversion; negative: femoral retroversion.
	- TT: positive: external rotation; negative: internal rotation.
note2: Adjust the MarkerSet.xml in the main folder to your marker set,
	when using markers for the greater trochanter (when adjusting
	femur) and/or when using markers on the feet (when adjusting tibia).
note3: If you only wish to adjust the femoral geometry (and not the tibial
	torsion), set the input to the tibial torsion to 0 degrees (=default
	tibial torsion in generic femur).
note4: Default angles of the generic OpenSim model geometry should be
  measured with the same method (e.g. Hernandez, ...) which you use for your 
  partipants to ensure consistency.
note5: in models with WrapObjects it is MANDATORY to check muscle moment
  arms for each movement before running muscle specific analysis (e.g.
  Static Optimization). 
  It can happen that muscles pull through WrapObjects, this has to be
  avoided. One possible solution to tackle this problem is to resize the 
  WrapObjects.
  You can use the "checkMuscleMomentArms" script to identfiy problems.


applyTibiaTorsionToJointOffset = 0 is the original method where torsion
is applied via translation and rotation axis and not via body coordinate
system rotation. This method is not applicable with Rajagopal model
because it does not have these elements...

12/26/2022
 changes Elias Wallnoefer:
 Femur torsion should now work with RajagopalModel
 FinalModel = "leftNSA*"

 ! Tibia torsion does not yet work - hard coded line-numbers of XML need to
 be fixed in tibia.m and tibia_locationInParent-rotation retested + adapted for both models


14/03/2023
 changes by Willi Koller
    Femur and Tibia Torsion should now work with all models, testet with gait2392, Hamner, Rajagopal, Lernagopal
    WrapObjects in the proximal part of the femur (part which is rotated) are not supported yet, you have to adjust the location and rotation manually.
    A message is written to the console if this is the case!

15/03/2023
 changes by Willi Koller for Lenhart model
    Lenhart model - need to adjust Ligament locations!
    also rotate additional geometries of foot

16/03/2023
 changes by Willi Koller for Lenhart model
    should work now with ligaments and additional geometries