/*Copyright 2022 Guillaume Villeneuve

   Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.*/

using Autodesk.AutoCAD.ApplicationServices;
using Autodesk.AutoCAD.Colors;
using Autodesk.AutoCAD.DatabaseServices;
using Autodesk.AutoCAD.EditorInput;
using Autodesk.AutoCAD.Geometry;
using Autodesk.AutoCAD.Runtime;
using System;
using System.Collections;
using System.IO;
using System.Windows.Forms;
using Application = Autodesk.AutoCAD.ApplicationServices.Application;


[assembly: ExtensionApplication(typeof(SACE_Traj_Gen.SACEPlugin))]

namespace SACE_Traj_Gen
{

    public class SACEPlugin : IExtensionApplication
    {
        // Plugin initialization
        #region Initialization
        void IExtensionApplication.Initialize()
        {
            Document doc = Autodesk.AutoCAD.ApplicationServices.Application.DocumentManager.MdiActiveDocument;
            Editor ed = doc.Editor;
            CreateSACELayer();
            ed.WriteMessage("\nSACE Plugin Loaded");

        }

        void IExtensionApplication.Terminate()
        {

        }
        #endregion
        // Useful subs
        #region Subs
        private static void CreateSACELayer(string saceLayerName = "SACE")
        {
            // Get the current document and database
            Document Doc = Application.DocumentManager.MdiActiveDocument;
            Database Db = Doc.Database;

            // Start a transaction
            using (Transaction trans = Db.TransactionManager.StartTransaction())
            {
                // Open the Layer table for read
                LayerTable acLyrTbl;
                acLyrTbl = trans.GetObject(Db.LayerTableId,
                                                OpenMode.ForRead) as LayerTable;

                if (acLyrTbl.Has(saceLayerName) == false)
                {
                    using (LayerTableRecord acLyrTblRec = new LayerTableRecord())
                    {
                        // Assign the layer the ACI color 3 and a name
                        acLyrTblRec.Color = Color.FromColorIndex(ColorMethod.ByAci, 3);
                        acLyrTblRec.Name = saceLayerName;

                        // Upgrade the Layer table for write
                        acLyrTbl.UpgradeOpen();

                        // Append the new layer to the Layer table and the transaction
                        acLyrTbl.Add(acLyrTblRec);
                        trans.AddNewlyCreatedDBObject(acLyrTblRec, true);
                    }
                }

                // Open the Block table for read
                BlockTable acBlkTbl;
                acBlkTbl = trans.GetObject(Db.BlockTableId,
                                                OpenMode.ForRead) as BlockTable;

                // Open the Block table record Model space for write
                BlockTableRecord acBlkTblRec;
                acBlkTblRec = trans.GetObject(acBlkTbl[BlockTableRecord.ModelSpace],
                                                OpenMode.ForWrite) as BlockTableRecord;

                // Save the changes and dispose of the transaction
                trans.Commit();
            }
        }

        #endregion
    }

    public static class SACETrajList
    // List of SACETraj objects in order
    {
        public static ArrayList SACETrajectories = new ArrayList();
        public static string folderpath = "";
        public static double max_error = 0.1;
        public static int traj_count = 0;
        public static double zlevel = 0.0;
        public static bool dynamic_point_density = true;
        public static double static_seg_len = 1;

        [CommandMethod("SACEToggleDynamicDensity")]
        public static void SetDynamicDensity()
        {
            Document doc = Autodesk.AutoCAD.ApplicationServices.Application.DocumentManager.MdiActiveDocument;
            Editor ed = doc.Editor;

            dynamic_point_density = !dynamic_point_density;

            ed.WriteMessage("Dynamic point density set to " + dynamic_point_density);
        }

        [CommandMethod("SACESetMaxError")]
        public static void SetMaxError()
        {
            Document doc = Autodesk.AutoCAD.ApplicationServices.Application.DocumentManager.MdiActiveDocument;
            Editor ed = doc.Editor;

            PromptDoubleResult maxerrorRes = ed.GetDouble("\nMaximum spline interpolation error: ");
            max_error = maxerrorRes.Value;
            ed.WriteMessage("Maximum spline interpolation value set to " + max_error);
        }

        [CommandMethod("SACESetStaticSegmentLength")]
        public static void SetStaticSegmentLength()
        {
            Document doc = Autodesk.AutoCAD.ApplicationServices.Application.DocumentManager.MdiActiveDocument;
            Editor ed = doc.Editor;

            PromptDoubleResult maxerrorRes = ed.GetDouble("\nStatic mode maximum segment length: ");
            static_seg_len = maxerrorRes.Value;
            ed.WriteMessage("Static mode maximum segment length " + static_seg_len);
        }

        [CommandMethod("SACESetMachiningPlane")]
        public static void SetMachiningPlance()
        {
            Document doc = Autodesk.AutoCAD.ApplicationServices.Application.DocumentManager.MdiActiveDocument;
            Editor ed = doc.Editor;

            PromptDoubleResult zlevelres = ed.GetDouble("\nMachining plane Z-level: ");
            zlevel = zlevelres.Value;
            ed.WriteMessage("Machining plane Z-level value set to " + zlevel);
        }



        [CommandMethod("SACEPathEngraving")]
        public static void PathToTrajectory()
        // Add a trajectory from linear geometry
        {
            Document doc = Autodesk.AutoCAD.ApplicationServices.Application.DocumentManager.MdiActiveDocument;
            Editor ed = doc.Editor;

            // TODO: Select all lines and join them
            // For now we merge manually

            PromptEntityOptions peo1 = new PromptEntityOptions("\nSelect trajectory geometry : ");
            peo1.SetRejectMessage("\nInvalid selection");
            PromptEntityResult pEntrs = ed.GetEntity(peo1);

            if (PromptStatus.OK != pEntrs.Status)
            {
                return;
            }

            // Selected object ID
            ObjectId lineID = pEntrs.ObjectId;

            // Create and add corresponding trajectory
            traj_count += 1;
            CurveToTraj(lineID, "traj" + traj_count + "_path");

        }

        [CommandMethod("SACEOffsetEngraving")]
        public static void OffsetToTrajectory()
        // Add a trajectory from linear geometry
        {
            Document doc = Autodesk.AutoCAD.ApplicationServices.Application.DocumentManager.MdiActiveDocument;
            Editor ed = doc.Editor;

            // TODO: Select all lines and join them
            // For now we merge manually

            PromptEntityOptions peo1 = new PromptEntityOptions("\nSelect trajectory geometry : ");
            peo1.SetRejectMessage("\nInvalid selection");
            PromptEntityResult pEntrs = ed.GetEntity(peo1);

            PromptDoubleResult resValueOff = ed.GetDistance("\nMaterial removal radius: ");
            // checks
            if (resValueOff.Status != PromptStatus.OK) return;
            PromptPointResult resPointOff = ed.GetPoint("\nOffset side: ");
            if (resPointOff.Status != PromptStatus.OK) return;

            // Create and add corresponding offset trajectory
            traj_count += 1;
            OffsetCurve(pEntrs, resValueOff.Value, resPointOff.Value, "traj" + traj_count + "_offset");

        }

        [CommandMethod("SACEConstantOffsets")]
        public static void ConstantOffsets()
        // Add a trajectory from linear geometry
        {
            Document doc = Autodesk.AutoCAD.ApplicationServices.Application.DocumentManager.MdiActiveDocument;
            Editor ed = doc.Editor;

            PromptEntityOptions peo1 = new PromptEntityOptions("\nSelect trajectory geometry : ");
            peo1.SetRejectMessage("\nInvalid selection");
            PromptEntityResult pEntrs = ed.GetEntity(peo1);

            PromptDoubleResult toolRadiusRes = ed.GetDistance("\nMaterial removal radius: ");
            if (toolRadiusRes.Status != PromptStatus.OK) return;

            PromptDoubleResult stepRes = ed.GetDistance("\nStep value: ");
            if (stepRes.Status != PromptStatus.OK) return;

            PromptDoubleResult iniRes = ed.GetDistance("\nInitial Offset: ");
            if (iniRes.Status != PromptStatus.OK) return;

            PromptPointResult resPointOff = ed.GetPoint("\nOffset side: ");
            if (resPointOff.Status != PromptStatus.OK) return;

            double radius = toolRadiusRes.Value;
            double step = stepRes.Value;
            double ini = iniRes.Value;

            traj_count += 1;
            int ofst_count = 1;

            for (double d = radius + ini; d > radius; d -= step)
            {
                double d2;

                if (d < radius)
                {
                    d2 = radius;
                }
                else
                {
                    d2 = d;
                }
                if (!(OffsetCurve(pEntrs, d2, resPointOff.Value, "traj" + traj_count + "_offset" + ofst_count)))
                {
                    return;
                }
                else ofst_count += 1;
            }
        }

        [CommandMethod("SACEDeepEngraving")]
        public static void DeepTrajectory()
        // Add a trajectory from linear geometry
        {
            Document doc = Autodesk.AutoCAD.ApplicationServices.Application.DocumentManager.MdiActiveDocument;
            Editor ed = doc.Editor;

            // TODO: Select all lines and join them
            // For now we merge manually

            PromptEntityOptions peo1 = new PromptEntityOptions("\nSelect trajectory geometry : ");
            peo1.SetRejectMessage("\nInvalid selection");
            PromptEntityResult pEntrs = ed.GetEntity(peo1);
            if (pEntrs.Status != PromptStatus.OK) return;

            PromptDoubleResult stepRes = ed.GetDouble("\nStep value: ");
            if (stepRes.Status != PromptStatus.OK) return;

            PromptDoubleResult final_depthRes = ed.GetDouble("\nFinal Depth: ");
            if (final_depthRes.Status != PromptStatus.OK) return;

            double step = stepRes.Value;
            double final_depth = final_depthRes.Value;

            // Selected object ID
            ObjectId lineID = pEntrs.ObjectId;

            // Create and add corresponding trajectory
            traj_count += 1;
            int dp_count = 1;

            for (double d = 0; d <= final_depth; d += step)
            {
                if (d <= final_depth)
                {
                    CurveToTraj(lineID, "traj" + traj_count + "_depth" + dp_count, d);
                    dp_count += 1;
                }
                else
                {
                    CurveToTraj(lineID, "traj" + traj_count + "_depth" + d, final_depth);
                    dp_count += 1;
                }

            }

        }


        //[CommandMethod("SACEHelicalDrill")]
        public static void HelicalDrill()
        // Add a trajectory to helical drill
        {

        }

        [CommandMethod("SACEDeleteLastTrajectory")]
        public static void DeleteTrajectory()
        // Remove a trajectory from the list
        {
            Document doc = Application.DocumentManager.MdiActiveDocument;
            Editor ed = doc.Editor;
            try
            {
                string trname = ((SACETraj)SACETrajectories[SACETrajectories.Count - 1]).name;
                SACETrajectories.RemoveAt(SACETrajectories.Count - 1);
                ed.WriteMessage("\n" + trname + " deleted");
            }
            catch
            {
                ed.WriteMessage("\nNo trajectory deleted");
            }
            ed.WriteMessage("\n" + SACETrajectories.Count + " trajectories remaining");
        }

        [CommandMethod("SACEDeleteAllTrajectories")]
        public static void DeleteAllTrajectories()
        // Remove a trajectory from the list
        {
            Document doc = Application.DocumentManager.MdiActiveDocument;
            Editor ed = doc.Editor;
            SACETrajectories.Clear();
            ed.WriteMessage("\n All trajectories cleared");
        }

        [CommandMethod("SACEExportTrajectories")]
        public static void ExportTrajectories()
        // Export all trajectories as .txt
        {
            Document doc = Application.DocumentManager.MdiActiveDocument;
            Editor ed = doc.Editor;
            FolderBrowserDialog folderBrowserDialog1 = new FolderBrowserDialog();
            if (folderBrowserDialog1.ShowDialog() == DialogResult.OK)
            {
                folderpath = folderBrowserDialog1.SelectedPath + "\\";

                foreach (SACETraj traj in SACETrajList.SACETrajectories)
                {
                    traj.ExportTrajectory();
                }
                ed.WriteMessage("\n{0} trajectories saved to {1}", SACETrajectories.Count, folderpath);
                SACETrajList.SACETrajectories.Clear();

            }
            else
            {
                ed.WriteMessage("export canceled");
            }
            
        }

        #region subs
        public static double GetLength(Curve ent)
        // Get length of Curve object
        {
            if (ent == null) return -1;

            return ent.GetDistanceAtParameter(ent.EndParam)
                - ent.GetDistanceAtParameter(ent.StartParam);
        }

        public static bool CurveToTraj(ObjectId lineID, string trajname = null, double depth = 0)
        // Save single curve object as SACETraj
        {
            Document doc = Application.DocumentManager.MdiActiveDocument;
            Editor ed = doc.Editor;
            Database db = doc.Database;

            // Start a transaction
            using (Transaction trans = db.TransactionManager.StartTransaction())
            {
                Curve TrajLine;
                try
                {
                    TrajLine = trans.GetObject(lineID, OpenMode.ForRead) as Curve;
                }
                catch
                {
                    ed.WriteMessage("\nInvalid selection");
                    return false;
                }

                if (TrajLine == null)
                {
                    ed.WriteMessage("\nInvalid selection");
                    return false;
                }

                SACETraj curTraj = new SACETraj();
                curTraj.length = GetLength(TrajLine);


                Point3d point0 = TrajLine.GetPointAtParameter(TrajLine.StartParam); // First point on curve
                Point3d point1 = TrajLine.GetPointAtParameter(TrajLine.EndParam); // Last point on curve
                Point3d pointdpi = TrajLine.GetPointAtDist(curTraj.length / 1E05); // point for initial derivative
                Point3d pointdpf = TrajLine.GetPointAtDist(curTraj.length - curTraj.length / 1E05); // point for final derivative

                double dxi = pointdpi.X - point0.X;
                double dyi = pointdpi.Y - point0.Y;
                double dzi = pointdpi.Z - point0.Z;
                double dxf = point1.X - pointdpf.X;
                double dyf = point1.Y - pointdpf.Y;
                double dzf = point1.Z - pointdpf.Z;

                // First point of spline to set ds
                Point3d pointini = new Point3d(point0.X - dxi, point0.Y - dyi, point0.Z - dzi);
                curTraj.pointList.Add(pointini);

                // First point on curve
                curTraj.pointList.Add(point0);

                //=Points Calculation===============================================================================
                double max_seg_len;
                double curLen = 0;
                Point3d newPt;

                double point_len = SACETrajList.max_error / 10;
                int point_num = (int)Math.Ceiling(curTraj.length / (point_len));



                Vector3d[] d1_list = new Vector3d[point_num];
                Vector3d[] d2_list = new Vector3d[point_num];
                double[] dd1_list = new double[point_num];
                double[] curv_list = new double[point_num];

                if (max_error > curTraj.length)
                {
                    return false;
                }

                if (dynamic_point_density)
                {
                    for (int i = 0; i < point_num; i++)
                    {
                        try
                        {
                            d1_list[i] = TrajLine.GetFirstDerivative(TrajLine.GetPointAtDist(i * point_len));
                            d2_list[i] = TrajLine.GetSecondDerivative(TrajLine.GetPointAtDist(i * point_len));
                            curv_list[i] = (d1_list[i].CrossProduct(d2_list[i]).Length) / Math.Pow(d1_list[i].Length, 3);
                        }
                        catch
                        {
                            d1_list[i] = d1_list[i - 1];
                            d2_list[i] = d2_list[i - 1];
                            curv_list[i] = 0;
                        }

                    };

                    for (int i = 1; i < point_num - 1; i++)
                    {
                        dd1_list[i] = (d1_list[i + 1].DivideBy(d1_list[i + 1].Length) - d1_list[i - 1].DivideBy(d1_list[i - 1].Length)).Length / (2 * point_len);
                    };
                }

                max_seg_len = curTraj.length / Math.Ceiling(curTraj.length / static_seg_len);

                while (curLen < curTraj.length)
                {
                    if (dynamic_point_density)
                    {
                        max_seg_len = GetMaxSegmentLength(TrajLine, curLen, dd1_list, curv_list, point_len, SACETrajList.max_error);
                    }

                    if (max_seg_len == 0) // A sharp corner (curvature = inf locally) was found
                    {
                        max_seg_len = SACETrajList.max_error / 8;

                        // Remove last point ; points are needed before and after 1st-order discontinuity
                        curTraj.pointList.RemoveAt(curTraj.pointList.Count - 1);
                        curLen -= max_error / 8;
                        if (curLen < 0)
                        {
                            curLen = 0;
                        }

                        for (int i = 0; i < 16; i++)
                        { //input 16 points to approximate corner spaced by max_error/8 ; due to properties of local control for cubic spline
                            if (curLen + max_seg_len < curTraj.length)
                            {
                                newPt = TrajLine.GetPointAtDist(curLen + max_seg_len);
                                curTraj.pointList.Add(newPt);
                                curLen += max_seg_len;
                            }
                            else
                            {
                                break;
                            }
                        }

                    }
                    else
                    {
                        if (curLen + max_seg_len < curTraj.length)
                        {
                            newPt = TrajLine.GetPointAtDist(curLen + max_seg_len);
                            curTraj.pointList.Add(newPt);
                            curLen += max_seg_len;
                        }
                        else
                        {
                            break;
                        }
                    }


                }




                //==================================================================================================


                curTraj.pointList.Add(point1);

                // Last point to set ds
                Point3d pointfin = new Point3d(point1.X + dxf, point1.Y + dyf, point1.Z + dzf);
                curTraj.pointList.Add(pointfin);
                #endregion
                if (trajname == null)
                {
                    curTraj.name = "traj" + (SACETrajectories.Count + 1).ToString();
                }
                else
                {
                    curTraj.name = trajname;
                }

                curTraj.translateZ(zlevel - depth);

                int result = SACETrajectories.Add(curTraj);
                ed.WriteMessage("\nTrajectory added : " + curTraj.name);

                // Dispose of the transaction
                trans.Commit();

                curTraj.PrintSpline();
                return true;
            }

        }

        public static double GetMaxSegmentLength(Curve TrajLine, double curLen, double[] dd1_list, double[] curv_list, double point_len, double maxerror)
        {
            int cur_pt = (int)Math.Ceiling(curLen / point_len);
            double try_len = GetLength(TrajLine) - curLen;
            if (curLen == 0)
            {
                try_len /= 2;
            }

            double max_len = maxerror;
            double error = maxerror + 1;
            double curv_max;

            bool corner_found = false;
            int try_pt;

            while (error > maxerror || corner_found)
            {

                if (try_len / point_len < 2.0)
                {
                    return 0;
                }

                corner_found = false;
                try_pt = (int)Math.Ceiling((curLen + try_len) / (point_len));
                curv_max = 0;

                for (int i = cur_pt; i < try_pt; i++)
                {

                    if (curv_list[i] > curv_max)
                    {
                        curv_max = curv_list[i];
                    }

                    if (dd1_list[i] > curv_list[i] * Math.Sqrt(2))
                    {
                        if (try_len / point_len < 2.0)
                            return 0; // signal that a sharp corner was found
                        else
                            corner_found = true;
                    }

                };


                max_len = try_len;
                error = (curv_max * Math.Pow(try_len, 2) / 24); // Error bound considering maximum second derivative
                try_len /= 2;
            }






            return max_len;
        }

        public static bool OffsetCurve(PromptEntityResult pEntrs, double offset_dist, Point3d sidePt, string trajname = null)
        {
            Document doc = Autodesk.AutoCAD.ApplicationServices.Application.DocumentManager.MdiActiveDocument;
            Editor ed = doc.Editor;

            using (Transaction tr = doc.TransactionManager.StartTransaction())
            {
                Curve curve = tr.GetObject(pEntrs.ObjectId, OpenMode.ForRead) as Curve;
                if (curve != null)
                {
                    BlockTableRecord btr = tr.GetObject(curve.BlockId, OpenMode.ForWrite) as BlockTableRecord;
                    if (btr != null)
                    {
                        Point3d pDir = (Point3d)(Application.GetSystemVariable("VIEWDIR"));
                        if (pDir != null)
                        {
                            Point3d pWCS = sidePt.TransformBy(ed.CurrentUserCoordinateSystem);
                            double offset = IsRightDirection(curve, pWCS, pDir.GetAsVector()) ? offset_dist : -offset_dist;
                            DBObjectCollection curvCols = curve.GetOffsetCurves(offset);

                            int sub_count = 1;

                            foreach (DBObject obj in curvCols)
                            {
                                Curve subCurv = obj as Curve;
                                if (subCurv != null)
                                {
                                    btr.AppendEntity(subCurv);
                                    //tr.AddNewlyCreatedDBObject(subCurv, true);

                                    // Selected object ID
                                    ObjectId lineID = subCurv.ObjectId;

                                    // Create and add corresponding trajectory
                                    CurveToTraj(lineID, trajname + "-" + sub_count);
                                    sub_count += 1;

                                }
                            }
                        }
                    }
                }
                tr.Commit();
            }
            if (PromptStatus.OK != pEntrs.Status)
            {
                return false;
            }
            else
            {
                return true;
            }
        }


        public static bool IsRightDirection(Curve pCurv, Point3d p, Vector3d vDir)
        // Detect side of point
        {
            Vector3d vNormal = Vector3d.ZAxis;
            if (pCurv.IsPlanar)
            {
                Plane plane = pCurv.GetPlane();
                vNormal = plane.Normal;
                p = p.Project(plane, vDir);
            }
            Point3d pNear = pCurv.GetClosestPointTo(p, true);
            Vector3d vSide = p - pNear;
            Vector3d vDeriv = pCurv.GetFirstDerivative(pNear);
            if (vNormal.CrossProduct(vDeriv).DotProduct(vSide) < 0.0)
                return true;
            else
                return false;
        }

    }

    public class SACETraj
    // One continuous trajectory to machine with given parameters
    {
        public ArrayList pointList;
        public double length;
        public string name = "trajectory";
        public double speed = 10; // mm/s
        public double voltage = 30; // Volts
        public string mode = "engraving";
        public ObjectId Splineid;
        public Point3d StartPt;

        public SACETraj()
        {
            pointList = new ArrayList { };
        }

        public void ExportTrajectory()
        {
            Document doc = Autodesk.AutoCAD.ApplicationServices.Application.DocumentManager.MdiActiveDocument;
            Editor ed = doc.Editor;
            string filename = name + ".txt";

            ArrayList pointListStr = new ArrayList();
            for (int i = 0; i < pointList.Count; i++)
            {
                pointListStr.Add(pointList[i].ToString());
            }


            String[] ptarray = (String[])pointListStr.ToArray(typeof(string));
            for (int i = 0; i < ptarray.Length; i++)
            {
                ptarray[i] = ptarray[i].Remove(0, 1);
                ptarray[i] = ptarray[i].Remove(ptarray[i].Length - 1, 1);
                ptarray[i] = ptarray[i].Insert(ptarray[i].IndexOf(',') + 1, " ");
                ptarray[i] = ptarray[i].Insert(ptarray[i].LastIndexOf(',') + 1, " ");
            }
            File.WriteAllLines(SACETrajList.folderpath + filename, ptarray);
        }

        public void PrintSpline()
        {
            if (pointList.Count == 0)
            {
                return;
            }
            Document doc = Autodesk.AutoCAD.ApplicationServices.Application.DocumentManager.MdiActiveDocument;
            Editor ed = doc.Editor;
            Database db = doc.Database;
            // Start a transaction
            using (Transaction acTrans = db.TransactionManager.StartTransaction())
            {
                // Open the Block table for read
                BlockTable acBlkTbl;
                acBlkTbl = acTrans.GetObject(db.BlockTableId,
                                             OpenMode.ForRead) as BlockTable;

                // Open the Block table record Model space for write
                BlockTableRecord acBlkTblRec;
                acBlkTblRec = acTrans.GetObject(acBlkTbl[BlockTableRecord.ModelSpace],
                                                OpenMode.ForWrite) as BlockTableRecord;



                Point3d[] ptarray = (Point3d[])pointList.ToArray(typeof(Point3d));
                Point3dCollection in_points = new Point3dCollection();

                for (int i = 1; i < ptarray.Length - 1; i++)
                {
                    in_points.Add(ptarray[i]);
                }

                Vector3d tan_ini = ptarray[1] - (Point3d)pointList[0];
                Vector3d tan_fin = ptarray[ptarray.Length - 1] - ptarray[ptarray.Length - 2];

                Spline acSpline = new Spline(in_points, tan_ini, tan_fin, 3, 0.0);

                acSpline.SetDatabaseDefaults();

                // Add the new object to the block table record and the transaction
                acBlkTblRec.AppendEntity(acSpline);
                acTrans.AddNewlyCreatedDBObject(acSpline, true);
                acSpline.Layer = "SACE";

                // Save the new line to the database
                acTrans.Commit();

                Splineid = acSpline.Id;
            }

        }

        public void translateZ(double zdist)
        {
            Vector3d translate = new Vector3d(0, 0, zdist);
            Point3d p;
            for (int i = 0; i < pointList.Count; i++)
            {
                p = (Point3d)pointList[i];
                pointList[i] = p.Add(translate);
            }
        }


    }


}

