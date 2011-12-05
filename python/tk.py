from Tkinter import *
import tkMessageBox
import tkFileDialog
import os

def StringVarListIndex(list,item):
    for i in range(len(list)):
        if list[i].get()==item.get():
            return i
    raise ValueError, "list.index(x): x not in list"

class HRA:
    def __init__(self):
        self.tk = Tk()
        self.tk.resizable(width=FALSE, height=FALSE)
        #self.tk.minsize(width=400, height=600)
        #self.tk.config(background="White")
        
        try:
            pass
            #self.tk.option_readfile("optionDB")
        except:
            pass
        
        #self.tk.geometry("800x600+0+0")
        self.page=None     

        #set mode
        self.editData=False
        self.editSurvey=False
        self.basicMode=True

        # create menubar
        menubar = Menu(self.tk)
        filemenu = Menu(menubar, tearoff=0)
        filemenu.add_command(label="New Basic", command=self.BasicMode)
        #filemenu.add_command(label="New Advanced", command=self.AdvancedMode)
        #filemenu.add_command(label="Open")
        #filemenu.add_command(label="Save")
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=self.tk.quit)
        menubar.add_cascade(label="File", menu=filemenu)

        helpmenu = Menu(menubar, tearoff=0)
        helpmenu.add_command(label="User Guide", command=self.OnUserGuide)
        menubar.add_cascade(label="Help", menu=helpmenu)

        # display the menu
        self.tk.config(menu=menubar)
        self.Home()
        
        self.tk.mainloop()

    def OnUserGuide(self):
        os.system("start C:\\InVEST_2.1.0\\InVEST_2.1.0_Documentation.pdf")

    def BasicMode(self):
        self.basicMode=True
        self.Home()

    def AdvancedMode(self):
        self.basicMode=False
        self.Home()

    def Home(self):
        """
        Displays the basic or advanced menu based on the user mode.
        """
        if self.basicMode:
            self.BasicMenu()
        else:
            self.AdvancedMenu()

            
    def NewData(self):
        """
        Initializes new variables for data.
        """
        self.habitats=StringVar()
        self.stressors=StringVar()
        self.GetHabitatStressors()
   

    def NewSurvey(self):
        """
        Initializes new variables for a survey.
        """
        self.dataClassNames=[]
        self.dataClasses=0
        
        self.stressorIntensityNames=[]
        self.stressorIntensities=0
        
        self.stressorManagementNames=[]
        self.stressorManagements=0
        
        self.areaChangeNames=[]
        self.areaChanges=0
        
        self.structureChangeNames=[]
        self.structureChanges=0
        
        self.disturbanceFrequencyNames=[]
        self.disturbanceFrequencies=0
        
        self.overlapTimeNames=[]
        self.overlapTimes=0
        
        self.mortalityRateNames=[]
        self.mortalityRates=0
        
        self.recruitmentPatternNames=[]
        self.recruitmentPatterns=0
        
        self.connectivityNames=[]
        self.connectivities=0
        
        self.recoveryTimeNames=[]
        self.recoveryTimes=0


    #example data
    def DefaultSurvey(self):
        """
        Sets the survey to the WCVI example.
        """
        dataClassNames=["unknown","limited","adequate","best"]      
        stressorIntensityNames=["unknown","low","medium","high"]
        stressorManagementNames=["unknown","not effective","somewhat effective","very effective"]
        areaChangeNames=["unknown","low (0-20%)","medium (20-50%)","high (50-100%)"]
        structureChangeNames=["unknown","low (0-20%)","medium (20-50%)","high (50-100%)"]
        disturbanceFrequencyNames=["unknown","daily to weekly","several times per year","annually or less often"]
        overlapTimeNames=["unknown","0-4 months","4-8 months","8-12 months"]
        mortalityRateNames=["unknown", "low (0-20%)", "moderate (20-50%)", "high (>=80%)"]
        recruitmentPatternNames=["unknown","annual or more often","every 1-2 years","every 2+ years"]
        connectivityNames=["unknown","low dispersal (<10km)","medium dispersal (10-100km)","high dispersal (100km+)"]
        recoveryTimeNames=["unknown","0-1 years","1-10 years","10+ years"]
        
        self.LoadSurvey(dataClassNames,
                        stressorIntensityNames,
                        stressorManagementNames,
                        areaChangeNames,
                        structureChangeNames,
                        disturbanceFrequencyNames,
                        overlapTimeNames,
                        mortalityRateNames,
                        recruitmentPatternNames,
                        connectivityNames,
                        recoveryTimeNames)

    def DefaultData(self):
        """
        Sets the data to the WCVI default.
        """
        self.DefaultSurvey()
        
        habitatNames=["kelp","eelgrass","hard bottom","soft bottom"]
        habitatQuality=[1,1,2,2]
                  
        stressorNames=["Finfish Aquaculture","Shellfish Aquaculture","Comm Salmon (Troll)","Rec Fishing"]
        stressorQuality=[1,1,2,2]
                  
        stressorIntensity=[3,2,3,2]
        stressorIntensityQuality=[0,0,0,0]

        stressorManagement=[3,2,3,2]
        stressorManagementQuality=[0,0,0,0]

        stressorDistance=[300,250,150,100]

        areaChange=[2,2,1,1,3,3,2,1,0,0,0,0,0,0,0,0]
        areaChangeQuality=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

        structureChange=[2,2,1,1,3,3,2,1,1,1,1,1,1,1,1,1]
        structureChangeQuality=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

        disturbanceFrequency=[2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
        disturbanceFrequencyQuality=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

        overlapTime=[3,3,2,2,3,3,2,2,3,3,2,2,3,3,2,2]
        overlapTimeQuality=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

        mortalityRate=[1,1,0,0]
        mortalityRateQuality=[0,0,0,0]

        recruitmentPattern=[2,1,0,0]
        recruitmentPatternQuality=[0,0,0,0]

        connectivity=[2,1,0,0]
        connectivityQuality=[0,0,0,0]

        recoveryTime=[1,1,3,1]
        recoveryTimeQuality=[0,0,0,0]


        self.LoadData(habitatNames,habitatQuality,
                      stressorNames,stressorQuality,
                      stressorIntensity,stressorIntensityQuality,
                      stressorManagement,stressorManagementQuality,
                      stressorDistance,
                      areaChange,areaChangeQuality,
                      structureChange,structureChangeQuality,
                      disturbanceFrequency,disturbanceFrequencyQuality,
                      overlapTime,overlapTimeQuality,
                      mortalityRate,mortalityRateQuality,
                      recruitmentPattern,recruitmentPatternQuality,
                      connectivity,connectivityQuality,
                      recoveryTime,recoveryTimeQuality)
    
    def WCVI(self):
        """
        Sets the survey and data to the WCVI example.
        """
        self.DefaultSurvey()
        self.DefaultData()

    #I/O
    def ReadFile(self,filename):
        """
        """
        infile=open(filename,"r")
        survey,habitat,stressor,habitatStressor=infile.read().strip().split("\n\n")
        infile.close()

        #survey
        survey=survey.split("\n")
        row=survey.pop(0).split(",")
        dataClasses=int(row.pop(0))
        dataClassNames=row

        row=survey.pop(0).split(",")
        stressorIntensities=int(row.pop(0))
        stressorIntensityNames=row

        row=survey.pop(0).split(",")
        stressorManagements=int(row.pop(0))
        stressorManagementNames=row

        row=survey.pop(0).split(",")
        areaChanges=int(row.pop(0))
        areaChangeNames=row

        row=survey.pop(0).split(",")
        structureChanges=int(row.pop(0))
        structureChangeNames=row

        row=survey.pop(0).split(",")
        disturbanceFrequencies=int(row.pop(0))
        disturbanceFrequencyNames=row

        row=survey.pop(0).split(",")
        overlapTimes=int(row.pop(0))
        overlapTimeNames=row
        
        row=survey.pop(0).split(",")
        mortalityRates=int(row.pop(0))
        mortalityRateNames=row

        row=survey.pop(0).split(",")
        recruitmentPatterns=int(row.pop(0))
        recruitmentPatternNames=row
        
        row=survey.pop(0).split(",")
        connectivities=int(row.pop(0))
        connectivityNames=row

        row=survey.pop(0).split(",")
        recoveryTimes=int(row.pop(0))
        recoveryTimeNames=row
        
        self.LoadSurvey(dataClassNames,
                        stressorIntensityNames,
                        stressorManagementNames,
                        areaChangeNames,
                        structureChangeNames,
                        disturbanceFrequencyNames,
                        overlapTimeNames,
                        mortalityRateNames,
                        recruitmentPatternNames,
                        connectivityNames,
                        recoveryTimeNames)


        #habitat data
        habitat=habitat.split("\n")
        habitat.pop(0) #skip column headers

        habitatNames=[]
        habitatQuality=[]
        mortalityRate=[]
        mortalityRateQuality=[]
        recruitmentPattern=[]
        recruitmentPatternQuality=[]
        connectivity=[]
        connectivityQuality=[]
        recoveryTime=[]
        recoveryTimeQuality=[]
        for row in habitat:
            row=row.split(",")
            habitatNames.append(row[1])
            habitatQuality.append(int(row[2]))
            mortalityRate.append(int(row[3]))
            mortalityRateQuality.append(int(row[4]))
            recruitmentPattern.append(int(row[5]))
            recruitmentPatternQuality.append(int(row[6]))
            connectivity.append(int(row[7]))
            connectivityQuality.append(int(row[8]))
            recoveryTime.append(int(row[9]))
            recoveryTimeQuality.append(int(row[10]))
        
        #stressor data
        stressor=stressor.split("\n")
        stressor.pop(0) #skip header row

        stressorNames=[]
        stressorQuality=[]       
        stressorIntensity=[]
        stressorIntensityQuality=[]
        stressorManagement=[]
        stressorManagementQuality=[]
        stressorDistance=[]
        for row in stressor:
            row=row.split(",")
            stressorNames.append(row[1])
            stressorQuality.append(int(row[2]))
            stressorIntensity.append(int(row[3]))
            stressorIntensityQuality.append(int(row[4]))
            stressorManagement.append(int(row[5]))
            stressorManagementQuality.append(int(row[6]))
            stressorDistance.append(float(row[7]))

        #habitat-stressor data
        habitatStressor=habitatStressor.split("\n")
        habitatStressor.pop(0) #skip header

        areaChange=[]
        areaChangeQuality=[]
        structureChange=[]
        structureChangeQuality=[]
        disturbanceFrequency=[]
        disturbanceFrequencyQuality=[]
        overlapTime=[]
        overlapTimeQuality=[]
        for row in habitatStressor:
            row=row.split(",")
            areaChange.append(int(row[4]))
            areaChangeQuality.append(int(row[5]))
            structureChange.append(int(row[6]))
            structureChangeQuality.append(int(row[7]))
            disturbanceFrequency.append(int(row[8]))
            disturbanceFrequencyQuality.append(int(row[9]))
            overlapTime.append(int(row[10]))
            overlapTimeQuality.append(int(row[11]))

        self.LoadData(habitatNames,habitatQuality,
                      stressorNames,stressorQuality,
                      stressorIntensity,stressorIntensityQuality,
                      stressorManagement,stressorManagementQuality,
                      stressorDistance,
                      areaChange,areaChangeQuality,
                      structureChange,structureChangeQuality,
                      disturbanceFrequency,disturbanceFrequencyQuality,
                      overlapTime,overlapTimeQuality,
                      mortalityRate,mortalityRateQuality,
                      recruitmentPattern,recruitmentPatternQuality,
                      connectivity,connectivityQuality,
                      recoveryTime,recoveryTimeQuality)            
            
    def WriteFile(self,filename):
        outfile=open(filename,"w")

        #survey        
        outfile.write(self.dataClasses.get()+","+",".join(map(lambda x: x.get(), self.dataClassNames)))
        outfile.write("\n"+self.stressorIntensities.get()+","+",".join(map(lambda x: x.get(), self.stressorIntensityNames)))
        outfile.write("\n"+self.stressorManagements.get()+","+",".join(map(lambda x: x.get(), self.stressorManagementNames)))
        outfile.write("\n"+self.areaChanges.get()+","+",".join(map(lambda x: x.get(), self.areaChangeNames)))
        outfile.write("\n"+self.structureChanges.get()+","+",".join(map(lambda x: x.get(), self.structureChangeNames)))
        outfile.write("\n"+self.disturbanceFrequencies.get()+","+",".join(map(lambda x: x.get(), self.disturbanceFrequencyNames)))
        outfile.write("\n"+self.overlapTimes.get()+","+",".join(map(lambda x: x.get(), self.overlapTimeNames)))
        outfile.write("\n"+self.mortalityRates.get()+","+",".join(map(lambda x: x.get(), self.mortalityRateNames)))
        outfile.write("\n"+self.recruitmentPatterns.get()+","+",".join(map(lambda x: x.get(), self.recruitmentPatternNames)))
        outfile.write("\n"+self.connectivities.get()+","+",".join(map(lambda x: x.get(), self.connectivityNames)))
        outfile.write("\n"+self.recoveryTimes.get()+","+",".join(map(lambda x: x.get(), self.recoveryTimeNames)))
        outfile.write("\n\n")

         #habitat data
        outfile.write("Habitat ID,Habitat Name,Habitat DQ,Mortality,Mortality DQ,Recruitment,Recruitment DQ,Connectivity,Connectivity DQ,Regeneration,Regeneration DQ")
        for i in range(int(self.habitats.get())):
            outfile.write("\n"+str(i+0))
            outfile.write(","+self.habitatNames[i].get())
            outfile.write(","+str(self.habitatQuality[i].get()+0))
            outfile.write(","+str(StringVarListIndex(self.mortalityRateNames,self.mortalityRate[i])))
            outfile.write(","+str(self.mortalityRateQuality[i].get()+0))
            outfile.write(","+str(StringVarListIndex(self.recruitmentPatternNames,self.recruitmentPattern[i])))
            outfile.write(","+str(self.recruitmentPatternQuality[i].get()+0))
            outfile.write(","+str(StringVarListIndex(self.connectivityNames,self.connectivity[i])))
            outfile.write(","+str(self.connectivityQuality[i].get()+0))
            outfile.write(","+str(StringVarListIndex(self.recoveryTimeNames,self.recoveryTime[i])))
            outfile.write(","+str(self.recoveryTimeQuality[i].get()+0))
        outfile.write("\n\n")
        
        #stressor data
        outfile.write("Stressor ID,Stressor Name,Stressor DQ,Intensity,Intensity DQ,Management,Managment DQ,Buffer")
        for i in range(int(self.stressors.get())):
            outfile.write("\n"+str(i+0))
            outfile.write(","+self.stressorNames[i].get())
            outfile.write(","+str(self.stressorQuality[i].get()+0))
            outfile.write(","+str(StringVarListIndex(self.stressorIntensityNames,self.stressorIntensity[i])))
            outfile.write(","+str(self.stressorIntensityQuality[i].get()+0))
            outfile.write(","+str(StringVarListIndex(self.stressorManagementNames,self.stressorManagement[i])))
            outfile.write(","+str(self.stressorManagementQuality[i].get()+0))
            outfile.write(","+str(self.stressorDistance[i].get()))
        outfile.write("\n\n")

        #habitat-stressor data
        outfile.write("Habitat ID,Habitat Name,Stressor ID,Stressor Name,Area Change,Area Change DQ,Structure Change,Structure Change DQ,Disturbance Frequency,Disturbance Frequency DQ,Temporal Overlap,Temporal Overlap DQ")
        for i in range(int(self.habitats.get())):
            for j in range(int(self.stressors.get())):
                outfile.write("\n"+str(i+0))
                outfile.write(","+self.habitatNames[i].get())
                outfile.write(","+str(j+0))
                outfile.write(","+self.stressorNames[j].get())
                outfile.write(","+str(StringVarListIndex(self.areaChangeNames,self.areaChange[(i*int(self.stressors.get()))+j])+0))
                outfile.write(","+str(self.areaChangeQuality[i*int(self.stressors.get())+j].get()+0))
                outfile.write(","+str(StringVarListIndex(self.structureChangeNames,self.structureChange[(i*int(self.stressors.get()))+j])+0))
                outfile.write(","+str(self.structureChangeQuality[i*int(self.stressors.get())+j].get()+0))
                outfile.write(","+str(StringVarListIndex(self.disturbanceFrequencyNames,self.disturbanceFrequency[(i*int(self.stressors.get()))+j])+0))
                outfile.write(","+str(self.disturbanceFrequencyQuality[i*int(self.stressors.get())+j].get()+0))
                outfile.write(","+str(StringVarListIndex(self.overlapTimeNames,self.overlapTime[(i*int(self.stressors.get()))+j])+0))
                outfile.write(","+str(self.overlapTimeQuality[i*int(self.stressors.get())+j].get()+0))
        outfile.close()
        

    def LoadValues(self,dataClassNames,
             stressorIntensityNames,
             stressorManagementNames,
             areaChangeNames,
             structureChangeNames,
             disturbanceFrequencyNames,
             overlapTimeNames,
             mortalityRateNames,
             recruitmentPatternNames,
             connectivityNames,
             recoveryTimeNames,
             habitatNames,habitatQuality,
             stressorNames,stressorQuality,
             stressorIntensity,stressorIntensityQuality,
             stressorManagement,stressorManagementQuality,
             stressorDistance,
             areaChange,areaChangeQuality,
             structureChange,structureChangeQuality,
             disturbanceFrequency,disturbanceFrequencyQuality,
             overlapTime,overlapTimeQuality,
             mortalityRate,mortalityRateQuality,
             recruitmentPattern,recruitmentPatternQuality,
             connectivity,connectivityQuality,
             recoveryTime,recoveryTimeQuality):
        """
        Sets the survey and data to the specified values.
        """
        self.LoadSurvey(dataClassNames,
                        stressorIntensityNames,stressorManagementNames,
                        areaChangeNames,structureChangeNames,
                        disturbanceFrequencyNames,overlapTimeNames,
                        mortalityRateNames,recruitmentPatternNames,
                        connectivityNames,recoveryTimeNames)
        self.LoadData(habitatNames,habitatQuality,
                      stressorNames,stressorQuality,
                      stressorIntensity,stressorIntensityQuality,
                      stressorManagement,stressorManagementQuality,
                      stressorDistance,
                      areaChange,areaChangeQuality,
                      structureChange,structureChangeQuality,
                      disturbanceFrequency,disturbanceFrequencyQuality,
                      overlapTime,overlapTimeQuality,
                      mortalityRate,mortalityRateQuality,
                      recruitmentPattern,recruitmentPatternQuality,
                      connectivity,connectivityQuality,
                      recoveryTime,recoveryTimeQuality)

    def LoadSurvey(self,dataClassNames,
                   stressorIntensityNames,
                   stressorManagementNames,
                   areaChangeNames,
                   structureChangeNames,
                   disturbanceFrequencyNames,
                   overlapTimeNames,
                   mortalityRateNames,
                   recruitmentPatternNames,
                   connectivityNames,
                   recoveryTimeNames):
        """
        Sets the survey to the specified values.
        """
        #survey        
        self.dataClasses=StringVar()
        self.dataClasses.set(str(len(dataClassNames)))
        self.dataClassNames=[]
        for i,name in enumerate(dataClassNames):
            self.dataClassNames.append(StringVar())
            self.dataClassNames[i].set(name)
            
        self.stressorIntensities=StringVar()
        self.stressorIntensities.set(str(len(stressorIntensityNames)))
        self.stressorIntensityNames=[]
        for i,name in enumerate(stressorIntensityNames):
            self.stressorIntensityNames.append(StringVar())
            self.stressorIntensityNames[i].set(name)

        self.stressorManagements=StringVar()
        self.stressorManagements.set(str(len(stressorManagementNames)))
        self.stressorManagementNames=[]
        for i,name in enumerate(stressorManagementNames):
            self.stressorManagementNames.append(StringVar())
            self.stressorManagementNames[i].set(name)
            
        self.areaChanges=StringVar()
        self.areaChanges.set(str(len(areaChangeNames)))
        self.areaChangeNames=[]
        for i,name in enumerate(areaChangeNames):
            self.areaChangeNames.append(StringVar())
            self.areaChangeNames[i].set(name)

        self.structureChanges=StringVar()
        self.structureChanges.set(str(len(structureChangeNames)))
        self.structureChangeNames=[]
        for i,name in enumerate(structureChangeNames):
            self.structureChangeNames.append(StringVar())
            self.structureChangeNames[i].set(name)

        self.disturbanceFrequencies=StringVar()
        self.disturbanceFrequencies.set(str(len(disturbanceFrequencyNames)))
        self.disturbanceFrequencyNames=[]
        for i,name in enumerate(disturbanceFrequencyNames):
            self.disturbanceFrequencyNames.append(StringVar())
            self.disturbanceFrequencyNames[i].set(name)

        self.overlapTimes=StringVar()
        self.overlapTimes.set(str(len(overlapTimeNames)))
        self.overlapTimeNames=[]
        for i,name in enumerate(overlapTimeNames):
            self.overlapTimeNames.append(StringVar())
            self.overlapTimeNames[i].set(name)

        self.mortalityRates=StringVar()
        self.mortalityRates.set(str(len(mortalityRateNames)))
        self.mortalityRateNames=[]
        for i,name in enumerate(mortalityRateNames):
            self.mortalityRateNames.append(StringVar())
            self.mortalityRateNames[i].set(name)
            
        self.recruitmentPatterns=StringVar()
        self.recruitmentPatterns.set(str(len(recruitmentPatternNames)))
        self.recruitmentPatternNames=[]
        for i,name in enumerate(recruitmentPatternNames):
            self.recruitmentPatternNames.append(StringVar())
            self.recruitmentPatternNames[i].set(name)

        self.connectivities=StringVar()
        self.connectivities.set(str(len(connectivityNames)))
        self.connectivityNames=[]
        for i,name in enumerate(connectivityNames):
            self.connectivityNames.append(StringVar())
            self.connectivityNames[i].set(name)
                                     
        self.recoveryTimes=StringVar()
        self.recoveryTimes.set(str(len(recoveryTimeNames)))
        self.recoveryTimeNames=[]
        for i,name in enumerate(recoveryTimeNames):
            self.recoveryTimeNames.append(StringVar())
            self.recoveryTimeNames[i].set(name)


    def LoadData(self,habitatNames,habitatQuality,
                 stressorNames,stressorQuality,
                 stressorIntensity,stressorIntensityQuality,
                 stressorManagement,stressorManagementQuality,
                 stressorDistance,
                 areaChange,areaChangeQuality,
                 structureChange,structureChangeQuality,
                 disturbanceFrequency,disturbanceFrequencyQuality,
                 overlapTime,overlapTimeQuality,
                 mortalityRate,mortalityRateQuality,
                 recruitmentPattern,recruitmentPatternQuality,
                 connectivity,connectivityQuality,
                 recoveryTime,recoveryTimeQuality):

        self.habitats=StringVar()
        self.habitats.set(str(len(habitatNames)))
        
        self.habitatNames=[]
        for i,habitat in enumerate(habitatNames):
            self.habitatNames.append(StringVar())
            self.habitatNames[i].set(habitat)
            
        self.habitatQuality=[]
        for i,quality in enumerate(habitatQuality):
            self.habitatQuality.append(IntVar())
            self.habitatQuality[i].set(quality)

        self.stressors=StringVar()
        self.stressors.set(str(len(stressorNames)))

        self.stressorNames=[]            
        for i,stressor in enumerate(stressorNames):
            self.stressorNames.append(StringVar())
            self.stressorNames[i].set(stressor)

        self.stressorQuality=[]
        for i,quality in enumerate(stressorQuality):
            self.stressorQuality.append(IntVar())
            self.stressorQuality[i].set(quality)
 
        self.stressorIntensity=[]
        for i,intensity in enumerate(stressorIntensity):
            self.stressorIntensity.append(StringVar())
            self.stressorIntensity[i].set(self.stressorIntensityNames[intensity].get())
        
        self.stressorIntensityQuality=[]
        for i,quality in enumerate(stressorIntensityQuality):
            self.stressorIntensityQuality.append(IntVar())
            self.stressorIntensityQuality[i].set(quality)
                                                
        self.stressorManagement=[]
        for i,management in enumerate(stressorManagement):
            self.stressorManagement.append(StringVar())
            self.stressorManagement[i].set(self.stressorManagementNames[management].get())
        
        self.stressorManagementQuality=[]
        for i,quality in enumerate(stressorManagementQuality):
            self.stressorManagementQuality.append(IntVar())
            self.stressorManagementQuality[i].set(quality)
    
        self.stressorDistance=[]
        for i,distance in enumerate(stressorDistance):
            self.stressorDistance.append(DoubleVar())
            self.stressorDistance[i].set(distance)
            
        self.areaChange=[]
        for i,area in enumerate(areaChange):
            self.areaChange.append(StringVar())
            self.areaChange[i].set(self.areaChangeNames[area].get())
                                   
        self.areaChangeQuality=[]
        for i,quality in enumerate(areaChangeQuality):
            self.areaChangeQuality.append(IntVar())
            self.areaChangeQuality[i].set(quality)
            
        self.structureChange=[]
        for i,change in enumerate(structureChange):
            self.structureChange.append(StringVar())
            self.structureChange[i].set(self.structureChangeNames[change].get())
            
        self.structureChangeQuality=[]
        for i,quality in enumerate(structureChangeQuality):
            self.structureChangeQuality.append(IntVar())
            self.structureChangeQuality[i].set(quality)
                                               
        self.disturbanceFrequency=[]
        for i,frequency in enumerate(disturbanceFrequency):
            self.disturbanceFrequency.append(StringVar())
            self.disturbanceFrequency[i].set(self.disturbanceFrequencyNames[frequency].get())
            
        self.disturbanceFrequencyQuality=[]
        for i,quality in enumerate(disturbanceFrequencyQuality):
            self.disturbanceFrequencyQuality.append(IntVar())
            self.disturbanceFrequencyQuality[i].set(quality)
            
        self.overlapTime=[]
        for i,t in enumerate(overlapTime):
            self.overlapTime.append(StringVar())
            self.overlapTime[i].set(self.overlapTimeNames[t].get())
            
        self.overlapTimeQuality=[]
        for i,quality in enumerate(overlapTimeQuality):
            self.overlapTimeQuality.append(IntVar())
            self.overlapTimeQuality[i].set(quality)
            
        self.mortalityRate=[]
        for i,mortality in enumerate(mortalityRate):
            self.mortalityRate.append(StringVar())
            self.mortalityRate[i].set(self.mortalityRateNames[mortality].get())
        
        self.mortalityRateQuality=[]
        for i,quality in enumerate(mortalityRateQuality):
            self.mortalityRateQuality.append(IntVar())
            self.mortalityRateQuality[i].set(quality)
        
        self.recruitmentPattern=[]
        for i,recruitment in enumerate(recruitmentPattern):
            self.recruitmentPattern.append(StringVar())
            self.recruitmentPattern[i].set(self.recruitmentPatternNames[recruitment].get())
        
        self.recruitmentPatternQuality=[]
        for i,quality in enumerate(recruitmentPatternQuality):
            self.recruitmentPatternQuality.append(IntVar())
            self.recruitmentPatternQuality[i].set(quality)
        
        self.connectivity=[]
        for i,c in enumerate(connectivity):
            self.connectivity.append(StringVar())
            self.connectivity[i].set(self.connectivityNames[c].get())
        
        self.connectivityQuality=[]
        for i,quality in enumerate(connectivityQuality):
            self.connectivityQuality.append(IntVar())
            self.connectivityQuality[i].set(quality)
        
        self.recoveryTime=[]
        for i,t in enumerate(recoveryTime):
            self.recoveryTime.append(StringVar())
            self.recoveryTime[i].set(self.recoveryTimeNames[t].get())

        self.recoveryTimeQuality=[]
        for i,quality in enumerate(recoveryTimeQuality):
            self.recoveryTimeQuality.append(IntVar())
            self.recoveryTimeQuality[i].set(quality)         
    


    #GUI
    def BasicMenu(self):
        if not self.page == None:
            self.page.pack_forget()

        self.page = Frame(self.tk)

        Label(self.page, text="Choose your survey?").grid(row=0,column=0,columnspan=4)
        Button(self.page, text="Default", command=self.EditSurvey, width=7).grid(row=1,column=1,sticky=W)
        Button(self.page, text="Import", command=self.Import, width=7).grid(row=1,column=2,sticky=E)

        self.page.pack()        
    
    def AdvancedMenu(self):
        if not self.page == None:
            self.page.pack_forget()

        self.page = Frame(self.tk)

        Label(self.page, text="Survey", width=7).grid(row=0,column=0)        
        Button(self.page, text="New", command=self.GetNewSurvey, width=7).grid(row=1,column=0)
        Button(self.page, text="Save", width=7).grid(row=2,column=0)
        Button(self.page, text="Open", width=7).grid(row=3,column=0)
        Button(self.page, text="Default", command=self.DefaultSurvey, width=7).grid(row=4,column=0)

        Label(self.page, text="Data", width=7).grid(row=0,column=1)
        Button(self.page, text="New", command=self.NewData, width=7).grid(row=1,column=1)
        Button(self.page, text="Import", command=self.Import, width=7).grid(row=2,column=1)
        Button(self.page, text="Export", width=7, command=self.Export).grid(row=3,column=1)
        Button(self.page, text="Example", command=self.DefaultData, width=7).grid(row=4,column=1)

        self.page.pack()

    def EditSurvey(self):
        if not self.page == None:
            self.page.pack_forget()

        self.page = Frame(self.tk)

        Label(self.page, text="Would you like to edit the survey labels?").grid(row=0,column=0,columnspan=4)
        Button(self.page, text="Yes", command=self.GetDataClassNames, width=3).grid(row=1,column=1,sticky=W)
        Button(self.page, text="No", command=self.GetNewData, width=3).grid(row=1,column=2,sticky=E)

        self.page.pack()            

    #survey GUI
    def GetNewSurvey(self):
        self.editSurvey=False
        self.dataClasses=StringVar()
        self.stressorIntensities=StringVar()
        self.stressorManagements=StringVar()
        self.areaChanges=StringVar()
        self.structureChanges=StringVar()
        self.disturbanceFrequencies=StringVar()
        self.overlapTimes=StringVar()
        self.mortalityRates=StringVar()
        self.recruitmentPatterns=StringVar()
        self.connectivities=StringVar()
        self.recoveryTimes=StringVar()
        self.GetSurvey()

    def GetSurvey(self):
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the number of possible values for each survey paramter.").grid(row=0, column=0, columnspan=6)

        Label(self.page, text="Data Quality").grid(row=1,column=1,sticky=E)
        Entry(self.page, width=2, textvariable=self.dataClasses).grid(row=1,column=2)

        Label(self.page, text="Intensity").grid(row=2,column=1,sticky=E)
        Entry(self.page, width=2, textvariable=self.stressorIntensities).grid(row=2,column=2)

        Label(self.page, text="Management").grid(row=2,column=3,sticky=E)
        Entry(self.page, width=2, textvariable=self.stressorManagements).grid(row=2,column=4)

        Label(self.page, text="Area Change").grid(row=3,column=1,sticky=E)
        Entry(self.page, width=2, textvariable=self.areaChanges).grid(row=3,column=2)

        Label(self.page, text="Structure Change").grid(row=3,column=3,sticky=E)
        Entry(self.page, width=2, textvariable=self.structureChanges).grid(row=3,column=4)

        Label(self.page, text="Disturbance Frequency").grid(row=4,column=1,sticky=E)
        Entry(self.page, width=2, textvariable=self.disturbanceFrequencies).grid(row=4,column=2)

        Label(self.page, text="Time Overlap").grid(row=4,column=3,sticky=E)
        Entry(self.page, width=2, textvariable=self.overlapTimes).grid(row=4,column=4)

        Label(self.page, text="Mortality Rate").grid(row=5,column=1,sticky=E)
        Entry(self.page, width=2, textvariable=self.mortalityRates).grid(row=5,column=2)

        Label(self.page, text="Recruitment").grid(row=5,column=3,sticky=E)
        Entry(self.page, width=2, textvariable=self.recruitmentPatterns).grid(row=5,column=4)

        Label(self.page, text="Connectivity").grid(row=6,column=1,sticky=E)
        Entry(self.page, width=2, textvariable=self.connectivities).grid(row=6,column=2)

        Label(self.page, text="Recovery Time").grid(row=6,column=3,sticky=E)
        Entry(self.page, width=2, textvariable=self.recoveryTimes).grid(row=6,column=4)

        #controls
        Button(self.page, text="Previous", command=self.Home, width=8).grid(row=7,column=0,sticky=W)
        Button(self.page, text="Next",command=self.GetSurveyValidate, width=8).grid(row=7,column=5,sticky=E)

        self.page.pack()

    def GetSurveyValidate(self):
        try:
            if self.editSurvey:
                pass
            else:
                self.dataClassNames=[]
                for i in range(int(self.dataClasses.get())):
                    self.dataClassNames.append(StringVar())

                self.stressorIntensityNames=[]
                for i in range(int(self.stressorIntensities.get())):
                    self.stressorIntensityNames.append(StringVar())
                    
                self.stressorManagementNames=[]
                for i in range(int(self.stressorManagements.get())):
                    self.stressorManagementNames.append(StringVar())
                    
                self.areaChangeNames=[]
                for i in range(int(self.areaChanges.get())):
                    self.areaChangeNames.append(StringVar())
                    
                self.structureChangeNames=[]
                for i in range(int(self.structureChanges.get())):
                    self.structureChangeNames.append(StringVar())
                    
                self.disturbanceFrequencyNames=[]
                for i in range(int(self.disturbanceFrequencies.get())):
                    self.disturbanceFrequencyNames.append(StringVar())

                self.overlapTimeNames=[]                    
                for i in range(int(self.overlapTimes.get())):
                    self.overlapTimeNames.append(StringVar())

                self.mortalityRateNames=[]                    
                for i in range(int(self.mortalityRates.get())):
                    self.mortalityRateNames.append(StringVar())

                self.recruitmentPatternNames=[]                    
                for i in range(int(self.recruitmentPatterns.get())):
                    self.recruitmentPatternNames.append(StringVar())

                self.connectivityNames=[]                    
                for i in range(int(self.connectivities.get())):
                    self.connectivityNames.append(StringVar())

                self.recoveryTimeNames=[]                    
                for i in range(int(self.recoveryTimes.get())):
                    self.recoveryTimeNames.append(StringVar())
                    
                self.GetDataClassNames()
        except ValueError:
            tkMessageBox.showerror("Value Error","All values must be integers.")
            self.GetSurvey()

    def GetDataClassNames(self):
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the data quality levels in increasing order.").grid(row=0, column=0, columnspan=4)

        for i in range(int(self.dataClasses.get())):
            Label(self.page, text=str(i+1)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.dataClassNames[i]).grid(row=i+1, column=2, sticky=W)

        #controls
        Button(self.page, text="Previous", command=self.Home, width=8).grid(row=int(self.dataClasses.get())+1,column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetIntensityNames, width=8).grid(row=int(self.dataClasses.get())+1,column=3,sticky=E)

        self.page.pack()
       


    def GetIntensityNames(self):
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the intensity levels in increasing order.").grid(row=0, column=0, columnspan=4)

        for i in range(int(self.stressorIntensities.get())):
            Label(self.page, text=str(i+1)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.stressorIntensityNames[i]).grid(row=i+1, column=2, sticky=W)
            
        #controls
        Button(self.page, text="Previous", command=self.GetDataClassNames, width=8).grid(row=int(self.stressorIntensities.get())+1,column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetManagementNames, width=8).grid(row=int(self.stressorIntensities.get())+1,column=3,sticky=E)

        self.page.pack()
        

    def GetManagementNames(self):
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the management levels in increasing order.").grid(row=0, column=0, columnspan=4)

        for i in range(int(self.stressorManagements.get())):
            Label(self.page, text=str(i+1)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.stressorManagementNames[i]).grid(row=i+1, column=2, sticky=W)
            
        #controls
        Button(self.page, text="Previous", command=self.GetIntensityNames, width=8).grid(row=int(self.stressorManagements.get())+1,column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetAreaChangeNames, width=8).grid(row=int(self.stressorManagements.get())+1,column=3,sticky=E)

        self.page.pack()


    def GetAreaChangeNames(self):
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the area change levels in increasing order.").grid(row=0, column=0, columnspan=4)

        for i in range(int(self.areaChanges.get())):
            Label(self.page, text=str(i+1)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.areaChangeNames[i]).grid(row=i+1, column=2, sticky=W)
            

        #controls
        Button(self.page, text="Previous", command=self.GetManagementNames, width=8).grid(row=int(self.areaChanges.get())+1,column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetStructureChangeNames, width=8).grid(row=int(self.areaChanges.get())+1,column=3,sticky=E)

        self.page.pack()


    def GetStructureChangeNames(self):
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the structure change levels in increasing order.").grid(row=0, column=0, columnspan=4)

        for i in range(int(self.structureChanges.get())):
            Label(self.page, text=str(i+1)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.structureChangeNames[i]).grid(row=i+1, column=2, sticky=W)
            

        #controls
        Button(self.page, text="Previous", command=self.GetAreaChangeNames, width=8).grid(row=int(self.structureChanges.get())+1,column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetDisturbanceFrequencyNames, width=8).grid(row=int(self.structureChanges.get())+1,column=3,sticky=E)

        self.page.pack()


    def GetDisturbanceFrequencyNames(self):
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the disturbance frequency levels in increasing order.").grid(row=0, column=0, columnspan=4)

        for i in range(int(self.disturbanceFrequencies.get())):
            Label(self.page, text=str(i+1)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.disturbanceFrequencyNames[i]).grid(row=i+1, column=2, sticky=W)
            

        #controls
        Button(self.page, text="Previous", command=self.GetStructureChangeNames, width=8).grid(row=int(self.disturbanceFrequencies.get())+1,column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetOverlapTimeNames, width=8).grid(row=int(self.disturbanceFrequencies.get())+1,column=3,sticky=E)

        self.page.pack()


    def GetOverlapTimeNames(self):
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the time overlap levels in increasing order.").grid(row=0, column=0, columnspan=4)

        for i in range(int(self.overlapTimes.get())):
            Label(self.page, text=str(i+1)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.overlapTimeNames[i]).grid(row=i+1, column=2, sticky=W)
            

        #controls
        Button(self.page, text="Previous", command=self.GetDisturbanceFrequencyNames, width=8).grid(row=int(self.overlapTimes.get())+1,column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetMortalityRateNames, width=8).grid(row=int(self.overlapTimes.get())+1,column=3,sticky=E)

        self.page.pack()


    def GetMortalityRateNames(self):
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the mortality rate levels in increasing order.").grid(row=0, column=0, columnspan=4)

        for i in range(int(self.mortalityRates.get())):
            Label(self.page, text=str(i+1)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.mortalityRateNames[i]).grid(row=i+1, column=2, sticky=W)
            

        #controls
        Button(self.page, text="Previous", command=self.GetOverlapTimeNames, width=8).grid(row=int(self.mortalityRates.get())+1,column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetRecruitmentPatternNames, width=8).grid(row=int(self.mortalityRates.get())+1,column=3,sticky=E)

        self.page.pack()


    def GetRecruitmentPatternNames(self):
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the recruitment pattern levels in increasing order.").grid(row=0, column=0, columnspan=4)

        for i in range(int(self.recruitmentPatterns.get())):
            Label(self.page, text=str(i+1)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.recruitmentPatternNames[i]).grid(row=i+1, column=2, sticky=W)
            

        #controls
        Button(self.page, text="Previous", command=self.GetMortalityRateNames, width=8).grid(row=int(self.recruitmentPatterns.get())+1,column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetConnectivityNames, width=8).grid(row=int(self.recruitmentPatterns.get())+1,column=3,sticky=E)

        self.page.pack()


    def GetConnectivityNames(self):
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the connectivity levels in increasing order.").grid(row=0, column=0, columnspan=4)

        for i in range(int(self.connectivities.get())):
            Label(self.page, text=str(i+1)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.connectivityNames[i]).grid(row=i+1, column=2, sticky=W)
            

        #controls
        Button(self.page, text="Previous", command=self.GetRecruitmentPatternNames, width=8).grid(row=int(self.connectivities.get())+1,column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetRecoveryTimeNames, width=8).grid(row=int(self.connectivities.get())+1,column=3,sticky=E)

        self.page.pack()


    def GetRecoveryTimeNames(self):
        self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the recovery time levels in increasing order.").grid(row=0, column=0, columnspan=6)

        for i in range(int(self.recoveryTimes.get())):
            Label(self.page, text=str(i+1)).grid(row=i+1, column=1, sticky=E)
            Entry(self.page, width=30, textvariable=self.recoveryTimeNames[i]).grid(row=i+1, column=2, sticky=W)
            

        #controls
        Button(self.page, text="Previous", command=self.GetConnectivityNames, width=8).grid(row=int(self.recoveryTimes.get())+1,column=0,sticky=W)
        Button(self.page, text="Done", command=self.Home, width=8).grid(row=int(self.recoveryTimes.get())+1,column=3,sticky=E)

        self.page.pack()



    #data GUI
    def GetNewData(self):
        """
        """
        self.editData=False
        self.DefaultSurvey()
        
        self.habitats=StringVar()
        self.stressors=StringVar()

        self.GetHabitatStressors()        

            
    def GetHabitatStressors(self):
        #allow starting GUI here
        if not self.page==None:
            self.page.pack_forget()
            
        self.page = Frame(self.tk)
        Label(self.page, text="Specify the number of habitats and stressors.").grid(row=0, column=0, columnspan=4)

        #habitats
        Label(self.page, text="Habitats").grid(row=1,column=1,sticky=E)
        Entry(self.page, textvariable=self.habitats, width=2).grid(row=1,column=2,sticky=W)

        #stressors
        Label(self.page, text="Stressors").grid(row=2,column=1,sticky=E)
        Entry(self.page, textvariable=self.stressors, width=2).grid(row=2,column=2,sticky=W)

        #controls
        Button(self.page, text="Previous", command=self.Home, width=8).grid(row=3,column=0,sticky=W)
        Button(self.page, text="Next",command=self.GetHabitatStressorsValidate, width=8).grid(row=3,column=3,sticky=E)
        self.page.pack()

    def GetHabitatStressorsValidate(self):
        """
        Validates habitat and stressor numbers. Initializes variables.
        """
        try:
            habitats=int(self.habitats.get())
            stressors=int(self.stressors.get())
            if self.editData:
                pass                
            else:
                habitatNames=[""]*habitats
                habitatQuality=[0]*habitats
                stressorNames=[""]*stressors
                stressorQuality=[0]*stressors
                stressorIntensity=[0]*stressors
                stressorIntensityQuality=[0]*stressors
                stressorManagement=[0]*stressors
                stressorManagementQuality=[0]*stressors
                stressorDistance=[0]*stressors
                areaChange=[0]*habitats*stressors
                areaChangeQuality=[0]*habitats*stressors
                structureChange=[0]*habitats*stressors
                structureChangeQuality=[0]*habitats*stressors
                disturbanceFrequency=[0]*habitats*stressors
                disturbanceFrequencyQuality=[0]*habitats*stressors
                overlapTime=[0]*habitats*stressors
                overlapTimeQuality=[0]*habitats*stressors
                mortalityRate=[0]*habitats
                mortalityRateQuality=[0]*habitats
                recruitmentPattern=[0]*habitats
                recruitmentPatternQuality=[0]*habitats
                connectivity=[0]*habitats
                connectivityQuality=[0]*habitats
                recoveryTime=[0]*habitats
                recoveryTimeQuality=[0]*habitats
                self.LoadData(habitatNames,habitatQuality,
                              stressorNames,stressorQuality,
                              stressorIntensity,stressorIntensityQuality,
                              stressorManagement,stressorManagementQuality,
                              stressorDistance,
                              areaChange,areaChangeQuality,
                              structureChange,structureChangeQuality,
                              disturbanceFrequency,disturbanceFrequencyQuality,
                              overlapTime,overlapTimeQuality,
                              mortalityRate,mortalityRateQuality,
                              recruitmentPattern,recruitmentPatternQuality,
                              connectivity,connectivityQuality,
                              recoveryTime,recoveryTimeQuality)
            self.GetHabitatNames()
        except ValueError:
            tkMessageBox.showerror("Value Error","All values must be integers.")
            self.GetHabitatStressors()

    def GetHabitatNames(self):
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify names for the habitats.").grid(row=0, column=0, columnspan=5+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Habitat Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Data Quality").grid(row=1,column=3, columnspan=int(self.dataClasses.get()))

        for i in range(int(self.habitats.get())):
            Label(self.page, text=str(i+1)).grid(row=2+i,column=1,sticky=E)
            Entry(self.page, width=20, textvariable=self.habitatNames[i]).grid(row=2+i,column=2)
            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.habitatQuality[i]).grid(row=2+i,column=3+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.habitats.get()),column=3+k,sticky=N)

        #controls
        Button(self.page, text="Previous", command=self.GetHabitatStressors, width=8).grid(row=3+int(self.habitats.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetStressorNames, width=8).grid(row=3+int(self.habitats.get()),column=4+int(self.dataClasses.get()),sticky=E)

        self.page.pack()


    def GetStressorNames(self):
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify names for the stressors.").grid(row=0, column=0, columnspan=5+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Stressor Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Data Quality").grid(row=1,column=3, columnspan=int(self.dataClasses.get()))
        
        for i in range(int(self.stressors.get())):
            Label(self.page, text=str(i+1)).grid(row=2+i,column=1,sticky=E)
            Entry(self.page, width=20, textvariable=self.stressorNames[i]).grid(row=2+i,column=2)
            
            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.stressorQuality[i]).grid(row=2+i,column=3+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.stressors.get()),column=3+k,sticky=N)

        #controls
        Button(self.page, text="Previous", command=self.GetHabitatNames, width=8).grid(row=3+int(self.stressors.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetIntensity, width=8).grid(row=3+int(self.stressors.get()),column=4+int(self.dataClasses.get()),sticky=E)

        self.page.pack()

    def GetIntensity(self):
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify stressor intensity.").grid(row=0, column=0, columnspan=6+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Stressor Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Intensity").grid(row=1,column=3)
        Label(self.page, text="Data Quality").grid(row=1,column=4, columnspan=int(self.dataClasses.get()))
        
        for i in range(int(self.stressors.get())):
            Label(self.page, text=str(i+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.stressorNames[i].get()).grid(row=2+i,column=2,sticky=W)
            apply(OptionMenu,[self.page, self.stressorIntensity[i]]+map(lambda n: n.get(),self.stressorIntensityNames)).grid(row=2+i,column=3,sticky=N+E+S+W)
                                                     
            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.stressorIntensityQuality[i]).grid(row=2+i,column=4+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.stressors.get()),column=4+k,sticky=N)

        #controls
        Button(self.page, text="Previous", command=self.GetStressorNames, width=8).grid(row=3+int(self.stressors.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetManagement, width=8).grid(row=3+int(self.stressors.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()

    def GetManagement(self):
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify stressor management.").grid(row=0, column=0, columnspan=4+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Stressor Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Management").grid(row=1,column=3)
        Label(self.page, text="Data Quality").grid(row=1,column=4, columnspan=int(self.dataClasses.get()))
        
        for i in range(int(self.stressors.get())):
            Label(self.page, text=str(i+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.stressorNames[i].get()).grid(row=2+i,column=2,sticky=W)
            apply(OptionMenu,[self.page, self.stressorManagement[i]]+map(lambda n: n.get(),self.stressorManagementNames)).grid(row=2+i,column=3,sticky=N+E+S+W)
                                                     
            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.stressorManagementQuality[i]).grid(row=2+i,column=4+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.stressors.get()),column=4+k,sticky=N)

        #controls
        Button(self.page, text="Previous", command=self.GetIntensity, width=8).grid(row=3+int(self.stressors.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetDistance, width=8).grid(row=3+int(self.stressors.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()

    def GetDistance(self):
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify stressor buffer distance.").grid(row=0, column=0, columnspan=5)

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Stressor Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Distance").grid(row=1,column=3)
        
        for i in range(int(self.stressors.get())):
            Label(self.page, text=str(i+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.stressorNames[i].get()).grid(row=2+i,column=2,sticky=W)
            Entry(self.page, width=20, textvariable=self.stressorDistance[i]).grid(row=2+i,column=3)
                                                     
        #controls
        Button(self.page, text="Previous", command=self.GetManagement, width=8).grid(row=3+int(self.stressors.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetAreaChange, width=8).grid(row=3+int(self.stressors.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()

    def GetAreaChange(self):
        self.habitatIndex=0
        self.GetHabitatAreaChange()

    def GetPreviousHabitatAreaChange(self):
        if self.habitatIndex == 0:
            self.GetDistance()
        else:
            self.habitatIndex=self.habitatIndex-1
            self.GetHabitatAreaChange()

    def GetNextHabitatAreaChange(self):
        if self.habitatIndex+1 < int(self.habitats.get()):
            self.habitatIndex=self.habitatIndex+1
            self.GetHabitatAreaChange()
        else:
            self.GetStructureChange()

    def GetHabitatAreaChange(self):
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify the habitat's change in stressor area.").grid(row=0, column=0, columnspan=6+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Habitat Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Id").grid(row=1,column=3,sticky=E)
        Label(self.page, text="Stressor Name").grid(row=1,column=4,sticky=W)
        Label(self.page, text="Change in Area").grid(row=1,column=5)
        Label(self.page, text="Data Quality").grid(row=1,column=6, columnspan=int(self.dataClasses.get()))

        for i in range(int(self.stressors.get())):
            Label(self.page, text=str(self.habitatIndex+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.habitatNames[self.habitatIndex].get()).grid(row=2+i,column=2,sticky=W)
            Label(self.page, text=str(i+1)).grid(row=2+i,column=3,sticky=E)
            Label(self.page, text=self.stressorNames[i].get()).grid(row=2+i,column=4,sticky=W)
            apply(OptionMenu,[self.page, self.areaChange[(self.habitatIndex*int(self.stressors.get()))+i]]+map(lambda n: n.get(),self.areaChangeNames)).grid(row=2+i,column=5,sticky=N+E+S+W)

            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.areaChangeQuality[(self.habitatIndex*int(self.stressors.get()))+i]).grid(row=2+i,column=6+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.stressors.get()),column=6+k,sticky=N)

        #controls
        Button(self.page, text="Previous", command=self.GetPreviousHabitatAreaChange, width=8).grid(row=3+int(self.stressors.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetNextHabitatAreaChange, width=8).grid(row=3+int(self.stressors.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()
            

    def GetStructureChange(self):
        self.habitatIndex=0
        self.GetHabitatStructureChange()

    def GetPreviousHabitatStructureChange(self):
        if self.habitatIndex == 0:
            self.habitatIndex=int(self.habitats.get())-1
            self.GetHabitatAreaChange()
        else:
            self.habitatIndex=self.habitatIndex-1
            self.GetHabitatStructureChange()

    def GetNextHabitatStructureChange(self):
        if self.habitatIndex+1 < int(self.habitats.get()):
            self.habitatIndex=self.habitatIndex+1
            self.GetHabitatStructureChange()
        else:
            self.GetDisturbanceFrequency()

    def GetHabitatStructureChange(self):
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify the habitat's change in stressor structure.").grid(row=0, column=0, columnspan=6+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Habitat Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Id").grid(row=1,column=3,sticky=E)
        Label(self.page, text="Stressor Name").grid(row=1,column=4,sticky=W)
        Label(self.page, text="Change in Structure").grid(row=1,column=5)
        Label(self.page, text="Data Quality").grid(row=1,column=6, columnspan=int(self.dataClasses.get()))

        for i in range(int(self.stressors.get())):
            Label(self.page, text=str(self.habitatIndex+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.habitatNames[self.habitatIndex].get()).grid(row=2+i,column=2,sticky=W)
            Label(self.page, text=str(i+1)).grid(row=2+i,column=3,sticky=E)
            Label(self.page, text=self.stressorNames[i].get()).grid(row=2+i,column=4,sticky=W)
            apply(OptionMenu,[self.page, self.structureChange[(self.habitatIndex*int(self.stressors.get()))+i]]+map(lambda n: n.get(),self.structureChangeNames)).grid(row=2+i,column=5,sticky=N+E+S+W)

            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.structureChangeQuality[(self.habitatIndex*int(self.stressors.get()))+i]).grid(row=2+i,column=6+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.stressors.get()),column=6+k,sticky=N)

        #controls
        Button(self.page, text="Previous", command=self.GetPreviousHabitatStructureChange, width=8).grid(row=3+int(self.stressors.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetNextHabitatStructureChange, width=8).grid(row=3+int(self.stressors.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()



    def GetDisturbanceFrequency(self):
        self.habitatIndex=0
        self.GetHabitatDisturbanceFrequency()

    def GetPreviousHabitatDisturbanceFrequency(self):
        if self.habitatIndex == 0:
            self.habitatIndex=int(self.habitats.get())-1
            self.GetHabitatStructureChange()
        else:
            self.habitatIndex=self.habitatIndex-1
            self.GetHabitatDisturbanceFrequency()
            
    def GetNextHabitatDisturbanceFrequency(self):
        if self.habitatIndex+1 < int(self.habitats.get()):
            self.habitatIndex=self.habitatIndex+1
            self.GetHabitatDisturbanceFrequency()
        else:
            self.GetTimeOverlap()

    def GetHabitatDisturbanceFrequency(self):
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify the habitat's stressor natural disturbance frequency.").grid(row=0, column=0, columnspan=6+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Habitat Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Id").grid(row=1,column=3,sticky=E)
        Label(self.page, text="Stressor Name").grid(row=1,column=4,sticky=W)
        Label(self.page, text="Disturbance Frequency").grid(row=1,column=5)
        Label(self.page, text="Data Quality").grid(row=1,column=6, columnspan=int(self.dataClasses.get()))

        for i in range(int(self.stressors.get())):
            Label(self.page, text=str(self.habitatIndex+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.habitatNames[self.habitatIndex].get()).grid(row=2+i,column=2,sticky=W)
            Label(self.page, text=str(i+1)).grid(row=2+i,column=3,sticky=E)
            Label(self.page, text=self.stressorNames[i].get()).grid(row=2+i,column=4,sticky=W)
            apply(OptionMenu,[self.page, self.disturbanceFrequency[(self.habitatIndex*int(self.stressors.get()))+i]]+map(lambda n: n.get(),self.disturbanceFrequencyNames)).grid(row=2+i,column=5,sticky=N+E+S+W)

            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.disturbanceFrequencyQuality[(self.habitatIndex*int(self.stressors.get()))+i]).grid(row=2+i,column=6+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.stressors.get()),column=6+k,sticky=N)

        #controls
        Button(self.page, text="Previous", command=self.GetPreviousHabitatDisturbanceFrequency, width=8).grid(row=3+int(self.stressors.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetNextHabitatDisturbanceFrequency, width=8).grid(row=3+int(self.stressors.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()



    def GetTimeOverlap(self):
        self.habitatIndex=0
        self.GetHabitatTimeOverlap()

    def GetPreviousHabitatTimeOverlap(self):
        if self.habitatIndex == 0:
            self.habitatIndex=int(self.habitats.get())-1
            self.GetHabitatDisturbanceFrequency()
        else:
            self.habitatIndex=self.habitatIndex-1
            self.GetHabitatTimeOverlap()

    def GetNextHabitatTimeOverlap(self):
        if self.habitatIndex+1 < int(self.habitats.get()):
            self.habitatIndex=self.habitatIndex+1
            self.GetHabitatTimeOverlap()
        else:
            self.GetMortalityRate()

    def GetHabitatTimeOverlap(self):
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify the habitat's stressor overlap time.").grid(row=0, column=0, columnspan=6+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Habitat Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Id").grid(row=1,column=3,sticky=E)
        Label(self.page, text="Stressor Name").grid(row=1,column=4,sticky=W)
        Label(self.page, text="Overlap Time").grid(row=1,column=5)
        Label(self.page, text="Data Quality").grid(row=1,column=6, columnspan=int(self.dataClasses.get()))

        for i in range(int(self.stressors.get())):
            Label(self.page, text=str(self.habitatIndex+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.habitatNames[self.habitatIndex].get()).grid(row=2+i,column=2,sticky=W)
            Label(self.page, text=str(i+1)).grid(row=2+i,column=3,sticky=E)
            Label(self.page, text=self.stressorNames[i].get()).grid(row=2+i,column=4,sticky=W)
            apply(OptionMenu,[self.page, self.overlapTime[(self.habitatIndex*int(self.stressors.get()))+i]]+map(lambda n: n.get(),self.overlapTimeNames)).grid(row=2+i,column=5,sticky=N+E+S+W)

            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.overlapTimeQuality[(self.habitatIndex*int(self.stressors.get()))+i]).grid(row=2+i,column=6+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.stressors.get()),column=6+k,sticky=N)

        #controls
        Button(self.page, text="Previous", command=self.GetPreviousHabitatTimeOverlap, width=8).grid(row=3+int(self.stressors.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetNextHabitatTimeOverlap, width=8).grid(row=3+int(self.stressors.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()


    def GetMortalityRate(self):
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify mortality rate.").grid(row=0, column=0, columnspan=6+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Habitat Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Mortality").grid(row=1,column=3)
        Label(self.page, text="Data Quality").grid(row=1,column=4, columnspan=int(self.dataClasses.get()))
        
        for i in range(int(self.habitats.get())):
            Label(self.page, text=str(i+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.habitatNames[i].get()).grid(row=2+i,column=2,sticky=W)
            apply(OptionMenu,[self.page, self.mortalityRate[i]]+map(lambda n: n.get(),self.mortalityRateNames)).grid(row=2+i,column=3,sticky=N+E+S+W)
                                                     
            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.mortalityRateQuality[i]).grid(row=2+i,column=4+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.habitats.get()),column=4+k,sticky=N)

        #controls
        Button(self.page, text="Previous", command=self.GetHabitatTimeOverlap, width=8).grid(row=3+int(self.habitats.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetRecruitmentPattern, width=8).grid(row=3+int(self.habitats.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()

    def GetRecruitmentPattern(self):
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify recruitment pattern.").grid(row=0, column=0, columnspan=6+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Habitat Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Pattern").grid(row=1,column=3)
        Label(self.page, text="Data Quality").grid(row=1,column=4, columnspan=int(self.dataClasses.get()))
        
        for i in range(int(self.habitats.get())):
            Label(self.page, text=str(i+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.habitatNames[i].get()).grid(row=2+i,column=2,sticky=W)
            apply(OptionMenu,[self.page, self.recruitmentPattern[i]]+map(lambda n: n.get(),self.recruitmentPatternNames)).grid(row=2+i,column=3,sticky=N+E+S+W)
                                                     
            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.recruitmentPatternQuality[i]).grid(row=2+i,column=4+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.habitats.get()),column=4+k,sticky=N)

       #controls
        Button(self.page, text="Previous", command=self.GetMortalityRate, width=8).grid(row=3+int(self.habitats.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetConnectivity, width=8).grid(row=3+int(self.habitats.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()

    def GetConnectivity(self):
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify connectivity.").grid(row=0, column=0, columnspan=6+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Habitat Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Connectivity").grid(row=1,column=3)
        Label(self.page, text="Data Quality").grid(row=1,column=4, columnspan=int(self.dataClasses.get()))
        
        for i in range(int(self.habitats.get())):
            Label(self.page, text=str(i+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.habitatNames[i].get()).grid(row=2+i,column=2,sticky=W)
            apply(OptionMenu,[self.page, self.connectivity[i]]+map(lambda n: n.get(),self.connectivityNames)).grid(row=2+i,column=3,sticky=N+E+S+W)
                                                     
            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.connectivityQuality[i]).grid(row=2+i,column=4+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.habitats.get()),column=4+k,sticky=N)

        #controls
        Button(self.page, text="Previous", command=self.GetRecruitmentPattern, width=8).grid(row=3+int(self.habitats.get()),column=0,sticky=W)
        Button(self.page, text="Next", command=self.GetRecoveryTime, width=8).grid(row=3+int(self.habitats.get()),column=6+int(self.dataClasses.get()),sticky=E)

        self.page.pack()

    def GetRecoveryTime(self):
        self.page.pack_forget()
        self.page = Frame(self.tk)

        Label(self.page, text="Specify recovery time.").grid(row=0, column=0, columnspan=6+int(self.dataClasses.get()))

        Label(self.page, text="Id").grid(row=1,column=1,sticky=E)
        Label(self.page, text="Habitat Name").grid(row=1,column=2,sticky=W)
        Label(self.page, text="Recovery").grid(row=1,column=3)
        Label(self.page, text="Data Quality").grid(row=1,column=4, columnspan=int(self.dataClasses.get()))
        
        for i in range(int(self.habitats.get())):
            Label(self.page, text=str(i+1)).grid(row=2+i,column=1,sticky=E)
            Label(self.page, text=self.habitatNames[i].get()).grid(row=2+i,column=2,sticky=W)
            apply(OptionMenu,[self.page, self.recoveryTime[i]]+map(lambda n: n.get(),self.recoveryTimeNames)).grid(row=2+i,column=3,sticky=N+E+S+W)
                                                     
            for j in range(int(self.dataClasses.get())):
                Radiobutton(self.page, text="", value=j, variable=self.recoveryTimeQuality[i]).grid(row=2+i,column=4+j)

        for k,dataClass in enumerate(self.dataClassNames):
            Label(self.page, text=dataClass.get(), wraplength=1).grid(row=2+int(self.habitats.get()),column=4+k,sticky=N)

        #controls
        Button(self.page, text="Previous", command=self.GetConnectivity, width=8).grid(row=3+int(self.habitats.get()),column=0,sticky=W)
        if self.basicMode:
            Button(self.page, text="Export", command=self.Export, width=8).grid(row=3+int(self.habitats.get()),column=6+int(self.dataClasses.get()),sticky=E)
        else:
            pass

        self.page.pack()

    def Export(self):
        filename=tkFileDialog.asksaveasfilename(filetypes=[("Comma Separated Values",".csv"),("All Files","*")])
        if filename==None:
            pass
        else:
            self.WriteFile(filename)

    def Import(self):
        filename=tkFileDialog.askopenfilename(filetypes=[("Comma Separated Values",".csv"),("All Files","*")])
        if filename=="":
            pass
        else:
            self.ReadFile(filename)
            self.EditSurvey()

if __name__ == "__main__":
    HRA()