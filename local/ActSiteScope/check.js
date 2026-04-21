
        var appData = {"pdb": "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N", "hotspots": [], "ligand": "", "ligand_format": "pdb", "mcsa": [], "consensus": [], "dx": "", "mode": "individual", "ligand_name": null, "show_stability": false, "show_base_residues": false, "vis_mcsa": false, "vis_mpbind": false, "vis_p2rank": false, "vis_fpocket": false, "vis_apbs": false, "theme": "Light"};
        var viewer = null;

        function updateStatus(msg) { console.log("3D:"+msg); if(msg.includes("READY")){ document.getElementById('sys-mode').innerText = msg.split(" ")[0]; } }
        function updateTelemetry(ac, gl, mode, ligName) { 
            document.getElementById('atom-stat').innerText = ac+" atoms"; 
            document.getElementById('sys-mode').innerText = "READY ("+mode+")"; 
            var ligElem = document.getElementById('ligand-stat');
            if(ligElem) {
                ligElem.innerText = (ligName && ligName.length > 0) ? "LIGAND: " + ligName : "LIGAND: NONE";
            }
        }
        function getResiIds(resStr) { if(!resStr) return []; var m=String(resStr).match(/\d+/g); return m?m.map(Number):[]; }

        function addSimpleLabel(v, text, center, bgColor, fontColor) {
            if(!center) return;
            v.addLabel(text, {position:{x:center[0],y:center[1],z:center[2]}, backgroundColor:bgColor||"#222", fontColor:fontColor||"white", fontSize:13, fontFamily:"sans-serif", backgroundOpacity:0.9, padding:4});
        }

        function addEPLabel(v, val, center, resInfo) {
            if(!center||val===undefined) return;
            var isZero = Math.abs(val) <= 0.001;
            var bgCol = isZero ? "#FFFFFF" : (val > 0 ? "#2563eb" : "#ef4444");
            var txtCol = isZero ? "black" : "white";
            var text = resInfo + " | EP " + val.toFixed(2) + " kT/e";
            v.addLabel(text, {position:{x:center[0]+1.5,y:center[1]+1.5,z:center[2]+1.5}, backgroundColor:bgCol, fontColor:txtCol, fontSize:14, backgroundOpacity:0.95});
        }

        var getAdaptiveColorFunc = function(atoms) {
            var maxB=0; for(var i=0;i<atoms.length;i++) { if(atoms[i].b>maxB) maxB=atoms[i].b; }
            var isExp = maxB < 65;
            return function(atom) {
                if(!appData.show_stability) return '#CCCCCC';
                var v=atom.b;
                if(isExp) { if(v<=15) return '#0053D6'; if(v<=30) return '#65CBF3'; if(v<=45) return '#FFDB13'; return '#FF7D45'; }
                else { if(v>=90) return '#0053D6'; if(v>=70) return '#65CBF3'; if(v>=50) return '#FFDB13'; return '#FF7D45'; }
            };
        };

        function initViewer() {
            setTimeout(function() {
                try {
                    if(typeof $3Dmol === 'undefined') { updateStatus("CRITICAL: 3Dmol.js not loaded."); return; }
                    var gl=null;
                    try { var c=document.createElement("canvas"); gl=c.getContext("webgl")||c.getContext("experimental-webgl"); } catch(e) {}
                    var glStatus = gl ? "WebGL: OK" : "WebGL: FAILED";
                    
                    var viewerEle = $("#viewer3d");
                    viewer = $3Dmol.createViewer(viewerEle, {backgroundColor:"#f8fafc"});
                    viewer.resize(); // 초기 크기 강제 계산

                    if(!appData.pdb || appData.pdb.length < 10) { updateStatus("Warning: No PDB Data"); return; }

                    // symmetries 에러 방지를 위해 멀티모델 옵션 검토
                    var m = viewer.addModel(appData.pdb, "pdb", {keepH: true});
                    var atomCount = m.selectedAtoms().length;
                var adaptiveColorFunc = getAdaptiveColorFunc(m.selectedAtoms());
                var baseStyle = {cartoon:{colorfunc:adaptiveColorFunc, opacity:0.85}};
                if(appData.show_base_residues) { baseStyle.stick = {colorfunc:adaptiveColorFunc, radius:0.05, opacity:0.4}; }
                viewer.setStyle({model:0}, baseStyle);

                if(appData.vis_mcsa && appData.mcsa && appData.mcsa.length>0) {
                    var m_ids = appData.mcsa.map(r=>parseInt(r.res_num)).filter(id=>!isNaN(id));
                    viewer.addStyle({model:0,resi:m_ids}, {stick:{color:"#000",radius:0.2},sphere:{color:"#000",radius:0.35}});
                    appData.mcsa.forEach(mc => {
                        var labelPos = mc.center;
                        if(!labelPos) {
                            var selAtoms = viewer.selectedAtoms({model:0, resi:[parseInt(mc.res_num)]});
                            if(selAtoms && selAtoms.length > 0) {
                                var sx=0, sy=0, sz=0;
                                selAtoms.forEach(a=>{sx+=a.x;sy+=a.y;sz+=a.z;});
                                labelPos = [sx/selAtoms.length, sy/selAtoms.length, sz/selAtoms.length];
                            }
                        }
                        if(labelPos && appData.mode !== "docking") addSimpleLabel(viewer, mc.res_name+" "+mc.res_num, labelPos, "#222", "white");
                    });
                }

                if((appData.mode === "individual" || appData.mode === "consensus") && appData.vis_apbs && appData.dx && appData.dx.length > 100) {
                    try {
                        console.log("[APBS] Received Payload data length: " + appData.dx.length);
                        var binaryStr = atob(appData.dx);
                        var bytes = new Uint8Array(binaryStr.length);
                        for(var i=0; i<binaryStr.length; i++) { bytes[i] = binaryStr.charCodeAt(i); }
                        var rawDxStr = pako.ungzip(bytes, {to: 'string'});
                        
                        var voldata = new $3Dmol.VolumeData(rawDxStr, "dx");
                        console.log("[APBS] VolumeData parsed successfully. Adding surface...");
                        viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity: 0.85, volscheme: new $3Dmol.Gradient.RWB(-5, 5), voldata: voldata}, {model: 0});
                        console.log("[APBS] Surface added.");
                    } catch(e) { console.error("[APBS Render Error]", e); }
                }

                try {
                    var adaptiveColorFunc = getAdaptiveColorFunc(m.selectedAtoms());
                    var baseStyle = {cartoon:{colorfunc:adaptiveColorFunc, opacity:0.85};
                    if(appData.show_base_residues) { baseStyle.stick = {colorfunc:adaptiveColorFunc, radius:0.05, opacity:0.4}; }
                    viewer.setStyle({model:0}, baseStyle);

                    if(appData.mode === "individual") {
                        if(appData.hotspots && appData.hotspots.length>0) {
                            appData.hotspots.forEach(h => {
                                try {
                                    var isMPB=h.label&&(h.label.includes("MPBind")||h.label.includes("AI-Based Binding")), isP2R=h.label&&h.label.includes("P2Rank"), isFP=h.label&&h.label.includes("fPocket");
                                    if(isMPB&&!appData.vis_mpbind) return;
                                    if(isP2R&&!appData.vis_p2rank) return;
                                    if(isFP&&!appData.vis_fpocket) return;
                                    var ids=getResiIds(h.residues), col=isP2R?"#FF5722":(isFP?"#2ECC71":"#00008F");
                                    viewer.addStyle({model:0,resi:ids}, {stick:{color:col,radius:0.2},sphere:{color:col,radius:0.35});
                                    if(h.center) viewer.addSphere({center:{x:h.center[0],y:h.center[1],z:h.center[2]},radius:5.0,color:col,opacity:0.1,wireframe:true});
                                    if(h.detailed_residues&&h.detailed_residues.length>0) { h.detailed_residues.forEach(dr=>{ if(dr.center) addSimpleLabel(viewer,dr.name+" "+dr.num,dr.center,col,"white"); }); }
                                } catch (e1) { console.warn("Individual label err:", e1); }
                            });
                        }
                    } else if(appData.mode === "consensus") {
                        if(appData.consensus&&appData.consensus.length>0) {
                            appData.consensus.forEach(h=>{ 
                                try {
                                    var ids=[parseInt(h.ResNo)]; 
                                    viewer.addStyle({model:0,resi:ids},{stick:{color:"#FFD700",radius:0.2},sphere:{color:"#FFD700",radius:0.35}); 
                                    if(appData.vis_apbs && h.ep_value !== undefined && h.center) {
                                        addEPLabel(viewer, h.ep_value, h.center, h.ResName+" "+h.ResNo); 
                                    } else if(h.center) {
                                        addSimpleLabel(viewer, h.ResName+" "+h.ResNo, h.center, "#FFD700", "#000");
                                    }
                                } catch(e2) { console.warn("Consensus label err:", e2); }
                            });
                        }
                    } else if(appData.mode === "docking") {
                        var labeledResIds = new Set();
                        if(appData.vis_apbs && appData.consensus && appData.consensus.length > 0) {
                            appData.consensus.forEach(h => {
                                try {
                                    var ids = [parseInt(h.ResNo)];
                                    labeledResIds.add(parseInt(h.ResNo));
                                    viewer.addStyle({model:0, resi:ids}, {stick:{color:"#FFD700", radius:0.2}, sphere:{color:"#FFD700", radius:0.35});
                                    if(h.ep_value !== undefined && h.center) {
                                        addEPLabel(viewer, h.ep_value, h.center, h.ResName + " " + h.ResNo);
                                    } else if(h.center) {
                                        addSimpleLabel(viewer, h.ResName + " " + h.ResNo, h.center, "#FFD700", "#000");
                                    }
                                } catch(e3) { console.warn("Docking consensus err:", e3); }
                            });
                        }
                        if(appData.vis_mcsa && appData.mcsa && appData.mcsa.length > 0) {
                            var m_ids = appData.mcsa.map(r => parseInt(r.res_num)).filter(id => !isNaN(id));
                            viewer.addStyle({model:0,resi:m_ids}, {stick:{color:"#000",radius:0.2},sphere:{color:"#000",radius:0.35});
                            appData.mcsa.forEach(mc => {
                                try {
                                    var resNum = parseInt(mc.res_num);
                                    if (labeledResIds.has(resNum)) return;
                                    
                                    var labelPos = mc.center;
                                    if(!labelPos) {
                                        var sa = viewer.selectedAtoms({model:0, resi:[resNum]});
                                        if(sa && sa.length > 0) {
                                            var sx=0, sy=0, sz=0; 
                                            sa.forEach(a => { sx+=a.x; sy+=a.y; sz+=a.z; });
                                            labelPos = [sx/sa.length, sy/sa.length, sz/sa.length];
                                        }
                                    }
                                    if(labelPos) {
                                        if (appData.vis_apbs && mc.ep_value !== undefined) {
                                            addEPLabel(viewer, mc.ep_value, labelPos, mc.res_name + " " + mc.res_num);
                                        } else {
                                            addSimpleLabel(viewer, mc.res_name + " " + mc.res_num, labelPos, "#222", "white");
                                        }
                                    }
                                } catch(e4) { console.warn("MCSA label err:", e4); }
                            });
                        }
                    }
                } catch(globalModeErr) { console.warn("Global Mode blocks err:", globalModeErr); }

                if(appData.ligand && appData.ligand.length>20) {
                    try {
                        var ligFmt = appData.ligand_format || "pdb";
                        var ligM = viewer.addModel(appData.ligand, ligFmt);
                        viewer.setStyle({model:1}, {stick:{color:"#800080",radius:0.22});
                        if(appData.ligand_name) {
                            var la=ligM.selectedAtoms(); if(la.length>0) { var cx=0,cy=0,cz=0; la.forEach(a=>{cx+=a.x;cy+=a.y;cz+=a.z;}); addSimpleLabel(viewer,appData.ligand_name,[cx/la.length,cy/la.length,cz/la.length],"#800080","white"); }
                        }
                    } catch(ligErr) { console.warn("Ligand load error:", ligErr); }
                }
                viewer.zoomTo({model:0});
                viewer.render();
                setTimeout(function() { viewer.resize(); viewer.center({model:0}); viewer.zoomTo({model:0}); viewer.render(); }, 300);
                setTimeout(function() { viewer.resize(); viewer.render(); updateTelemetry(atomCount, glStatus, appData.mode, appData.ligand_name); $(".loading").fadeOut(); }, 1200);

                } catch(e) { console.error("Init Error:", e); updateStatus("ENGINE ERROR: "+e.message); }
            }, 150);
        }
        $(document).ready(function() { setTimeout(initViewer, 500); });
    