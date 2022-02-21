# SACE_Traj_Gen

An AutoCAD plugin to facilitate generation of spline control points from drawing geometry. Useful for CAM where G-Code cannot be used and spline trajectories are needed.

# Commands

SACEPATHENGRAVING
Generates a spline from a selected curve and adds it to the trajectory list

SACEOFFSETENGRAVING
Generates a spline with an offset defined by the user

SACECONSTANTOFFSETS
Generates consecutive evenly-spaced offset splines, with the last one being on tool-radius distance from the curve

SACEDEEPENGRAVING
Generates a spline with a defined Z depth, starting at the current machining plane

SACEEXPORTTRAJECTORIES
Saves all trajectories as .txt the clears trajectory list

SACEDELETELASTTRAJECTORY
Deletes last generated trajectory

SACEDELETEALLTRAJECTORIES
Deletes all trajectories in list

SACESETMACHININGINGPLANE
Sets the Z-level of generated splines

SACETOGGLEDYNAMICDENSITY
Toggles whether splines will have uniformly-spaced contorl points of if their density is variable

SACESETMAXERROR
Sets maximum acceptable error; used in calculating control point density

SACESETSTATICSEGMENTLENGTH
Sets distance between control points when dynamic density is set to false