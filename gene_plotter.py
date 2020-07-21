# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 17:10:16 2020

@author: Bastian Hornung
email bastian dot hornung at gmx dot germany-tld

Written in Python3

This script takes one or multiple genbank file, and makes a basic plot with pre-specified
genes out of it.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
import os

VERSION = "0.9"

class entry:
    def __init__(self,sType,sStart,sStop,bComp=False):
        self.sType = sType
        self.iStart = int(sStart)
        self.iStop = int(sStop)
        self.sProduct = ""
        self.sGeneName = ""
        self.sLocus = ""
        self.sOldLocus = ""
        self.sProtein = ""
        self.bCompl = bComp
        self.sColour = "grey"
        self.lExons = []
        self.lIntrons = []


def get_colour(newEntry,dicColor):
    """"very simple colouring scheme. Can be manually overriden by a given file with hex colours
    No error handling for wrong colours in the file
    """
    
    sAnno = newEntry.sProduct.lower()
    if newEntry.sLocus in dicColor:
        return dicColor.get(newEntry.sLocus)
    if newEntry.sOldLocus in dicColor:
        return dicColor.get(newEntry.sOldLocus)    
    if newEntry.sProtein in dicColor:
        return dicColor.get(newEntry.sProtein)      
    if newEntry.sGeneName in dicColor:
        return dicColor.get(newEntry.sGeneName)
    if newEntry.sProduct in dicColor:
        return dicColor.get(newEntry.sProduct)
    if newEntry.sType=="tRNA" or newEntry.sType=="rRNA":
        return "#FF00FF"
    if "transpos" in sAnno or "phage" in sAnno or "integras" in sAnno:
        return "red"
    if "flagel" in sAnno:
        return "#800000"
    if "transport" in sAnno or "abc " in sAnno or "pts " in sAnno or "mfs" in sAnno or "membrane" in sAnno:
        return "#008000" # green
    if "ase " in sAnno:
        return "blue"
    if "ribosom" in sAnno:
        return "#00FF00" # lime
    return "grey"

def hasNumbers(inputString):
    return any(char.isdigit() for char in inputString)

def sanitize_organism_name(sName):
    ''' organism names should be in italics, but strain identifiers etc should be not
    This function is supposed to take care of a few obvious things, but will probably not
    catch most cases, since it's a bit hard to say which item is a strain identifier or not

    Other note:
    #in italics, you need to escape space signs with \, else they do not show up. Not relevant for
    this function anymore, but noted in case it appears again
    '''
    
    sName = sName.split("=",1)[1].strip(chr(34))    
    lName = sName.split(" ")
    sNameNew = ""
    lNot = ["sp.","subsp."]
    lEnd = ["DSM","ATCC","str."]
    bEnd = False
    for item in lName:
        if item in lEnd:
            bEnd = True
        #logic: Items in the lNot list are not italics. Items after str., DSM or ATCC are not italics, including these items
        #items which contain a number are not italics. Because no name has numbes in them, right?
        #items which are only uppercase are not in italics either, because that has to be a strain identifier
        if bEnd or item in lNot or hasNumbers(item) or item==item.upper():
            sNameNew = sNameNew+" "+item            
        else:
            sNameNew = sNameNew+" "+'$\it{'+item+'}$'

    return sNameNew

def process_location(sType,bComp,sLocation):
    sLocation = sLocation.replace(">","")
    sLocation = sLocation.replace("<","")    
    sLocation = sLocation.strip(")")
    if not ".." in sLocation:
        sLocation = sLocation+".."+sLocation #this is just to read in a SNP, which can be found as "variation" and a single location only
    lCoord = sLocation.split("..")
    sStart = lCoord[0]
    sStop = lCoord[-1]
    cNewEntry = entry(sType,sStart,sStop,bComp)

    if len(lCoord)>2:
        lParts = sLocation.split(",")
        for parts in lParts:
            s1,s2 = parts.split("..")
            cNewEntry.lExons.append([int(s1),int(s2)])
        for i,exons in enumerate(cNewEntry.lExons):
            if i==len(cNewEntry.lExons)-1:break
            cNewEntry.lIntrons.append([cNewEntry.lExons[i][1],cNewEntry.lExons[i+1][0]])
    return cNewEntry

def read_genbank(sFile,dicColor):
    """Basic function to read genbank files.
    Does not parse everything, only what is relevant.
    I did not want to use bio-python, because in this way I am independent,
    and this function is not that complicated to write.
    Although it begins to get long now, but most things should be covered
    """
    print ("processing: ",sFile)
    lEntry = []
    inputFile = open(sFile)
    cNewEntry = ""
    bRead = False
    sName = ""
    for lines in inputFile:
        lines = lines.strip("\n")
        if lines.startswith("                     /organism="):
            sName = sanitize_organism_name(lines)
        if lines.startswith("                     /strain="):
            sStrain = lines.strip().split("=")[1].strip(chr(34))
            if not sStrain in sName:
                sName = sName+" str. "+sStrain
        if lines.startswith("                     /sub_strain="):
            sSub = lines.strip().split("=")[1].strip(chr(34))
            if not sSub in sName:
                sName = sName+" substr. "+sSub                
        if lines.startswith("     gene            "):
            bRead = True
        if not bRead:continue
        if lines.startswith("ORIGIN"):
            bRead = False
        if lines.startswith("     ") and not lines.startswith("                     ") and not "source" in lines:
            if cNewEntry:
                cNewEntry.sColour = get_colour(cNewEntry,dicColor)
                lEntry.append(cNewEntry)
            if lines.count("..")<=1 or lines.endswith(")"):
                bSingleLine = True
            else:
                bSingleLine = False
            lData = lines.split(" ")
            for item in lData:
                if item:
                    sType = item
                    break
            sLocation = lData[-1]
            bComp = False
            if "order" in sLocation:
                sLocation = sLocation.replace("order(","")
            if "join" in sLocation:
                sLocation = sLocation.replace("join(","")
                if "complement" in sLocation:
                    sLocation = sLocation.replace("))",")")
                else:
                    if not sLocation.endswith(","):
                        sLocation = sLocation[:-1]
            if "complement" in sLocation:
                sLocation = sLocation.split("(")[1]
                sLocation = sLocation.strip(")")
                bComp = True
            if bSingleLine:
                cNewEntry = process_location(sType,bComp,sLocation)
            else:
                continue

        if not bSingleLine:
            if lines.strip().startswith("/"):
                bSingleLine = True
                cNewEntry = process_location(sType,bComp,sLocation)
            else:
                sLocation = sLocation+lines.strip()
                
        lines = lines.strip()
        lines = lines.replace(chr(34),"")
        #I am accutely aware that this could be smarter
        if lines.startswith("/gene="):
            sGene = lines.split("=")[1]
            cNewEntry.sGeneName = sGene
        if lines.startswith("/product="):
            sProduct = lines.split("=")[1]
            cNewEntry.sProduct = sProduct
        if lines.startswith("/locus_tag="):
            sLocus = lines.split("=")[1]
            cNewEntry.sLocus = sLocus
        if lines.startswith("/old_locus_tag="):
            sOldLocus = lines.split("=")[1]
            cNewEntry.sOldLocus = sOldLocus                
        if lines.startswith("/protein_id="):
            sProtein = lines.split("=")[1]
            cNewEntry.sProtein = sProtein


    cNewEntry.sColour = get_colour(cNewEntry,dicColor)
    lEntry.append(cNewEntry)
    return lEntry,sName

def get_label(sLabel,item,dicNames):
    '''creates the final output string. If anything is specified in the name file, it will be used.
    Otherwise the type of label can be specified, and the name will be composed. Default is locus+product
    Gene names will be in italics.
    '''
    
    if item.sLocus in dicNames:
        return dicNames.get(item.sLocus)
    if item.sOldLocus in dicNames:
        return dicNames.get(item.sOldLocus)    
    if item.sProtein in dicNames:
        return dicNames.get(item.sProtein)       
    if item.sGeneName in dicNames:
        return dicNames.get(item.sGeneName)
    if item.sProduct in dicNames:
        return dicNames.get(item.sProduct)
    if sLabel=="gene_name" and item.sGeneName:
        return '$\it{'+item.sGeneName+'}$'
    if sLabel=="locus" and item.sLocus:
        return item.sLocus
    if sLabel=="product" and item.sProduct:
        return item.sProduct
    if sLabel=="locus+gene_name" and item.sGeneName: #this is intentional, also below
        return item.sLocus+" "+'$\it{'+item.sGeneName+'}$'
    if sLabel=="gene_name+product" and item.sGeneName:
        return '$\it{'+item.sGeneName+'}$'+" "+item.sProduct
    if sLabel=="locus+product" and (item.sLocus or item.sProduct):
        return item.sLocus+" "+item.sProduct
    return item.sLocus#+" "+item.sProduct

def make_plot(lEntry,sRev,iScale,sLabel,sLabelPos,iRotation,sEntryType,sStartGene,sStopGene,sOut,sExt,iStartCoord,iStopCoord,iDistOffset,iSizeText,dicNames,sOrgName,ax,bCoord,fThick):
    '''plots each entry in the input file
    I should consider splitting up this function, it is a bit long and complicated
    '''
    
    iStopCoordOrg = iStopCoord
    iStopCoord = iStartCoord+iScale
    plt.xlim(iStartCoord-200,iStopCoord+200)
    
    #I admit, this below is maybe not the best understandable
    if sRev=="reverse":
        if bCoord:
            plt.annotate(iStopCoordOrg+200,(iStartCoord-200,0),fontsize=iSizeText,horizontalalignment="right")
            plt.annotate(sOrgName,(iStartCoord-200,0.04),fontsize=iSizeText,horizontalalignment="left")
            plt.annotate(iStartCoord-200,(iStopCoordOrg+200,0),fontsize=iSizeText)
        iEndH = float(iStopCoordOrg-iStartCoord)/float(iScale)
    else: #this is more straight forward
        if bCoord:
            plt.annotate(iStartCoord-200,(iStartCoord-200,0),fontsize=iSizeText,horizontalalignment="right")
            plt.annotate(iStopCoordOrg+200,(iStopCoordOrg+200,0),fontsize=iSizeText)
            plt.annotate(sOrgName,(iStartCoord-200,0.04),fontsize=iSizeText,horizontalalignment="left")
        iEndH = float(iStopCoordOrg-iStartCoord+200)/float(iStopCoord-iStartCoord+200)
        
    plt.axhline(0,color="black",linewidth=2,xmax=iEndH)        
    plt.ylim(-0.1,0.1)
    plt.axis('off')
    fThickFin = 0.3 * fThick
    #arrow_style="simple,head_length=0.1,head_width=0.3,tail_width=0.3"
    arrow_style="simple,head_length=0.1,head_width="+str(fThickFin)+",tail_width="+str(fThickFin)

    bPrint = False
    for item in lEntry:
        if (item.sLocus ==sStartGene or item.sGeneName==sStartGene or item.sProduct==sStartGene or item.sOldLocus==sStartGene or item.sProtein==sStartGene)  and (item.sType in sEntryType or not sEntryType):
            bPrint = True
        if bPrint and (item.sType in sEntryType or not sEntryType):
            if item.bCompl:
                iStartPrint = item.iStop
                iStopPrint = item.iStart
            else:
                iStartPrint = item.iStart
                iStopPrint = item.iStop
            arrow = mpatches.FancyArrowPatch((iStartPrint, 0), (iStopPrint, 0),mutation_scale=100,facecolor=item.sColour,arrowstyle=arrow_style)
            #arrow = mpatches.FancyArrowPatch((iStartPrint, 0), (iStartPrint+(iLen*9), 0),mutation_scale=100,facecolor=item.sColour,arrowstyle=arrow_style)
            ax.add_patch(arrow)
            iLength = max(item.iStop,item.iStart)-min(item.iStop,item.iStart)
            iMiddle = min(item.iStart,item.iStop)+iLength/2
            iY = 0.02
            if sLabelPos=="Down":iY = iY*-1
            iY = iY+iDistOffset
            sLabelOut = get_label(sLabel,item,dicNames)
            if not iRotation:
                custAlign="center"
            else:
                custAlign="left"
            plt.annotate(sLabelOut,(iMiddle,iY),rotation=iRotation,fontsize=iSizeText,horizontalalignment=custAlign)
            print (fThick,fThickFin,fig.get_size_inches()*fig.dpi)
            if item.lIntrons:                
                for introns in item.lIntrons:                    
                    if introns[1]-introns[0]<=0:continue #can't draw introns of negative size
                    rec = mpatches.Rectangle((introns[0],-0.019*fThick),width=introns[1]-introns[0],height=0.038*fThick,color="#A9A9A9")
                    ax.add_patch(rec)

        if (item.sLocus ==sStopGene or item.sGeneName==sStopGene or item.sProduct==sStopGene or item.sProtein==sStopGene or item.sOldLocus==sStopGene)  and (item.sType in sEntryType or not sEntryType):
            bPrint = False
            break

    return True

def get_start_stop_coords(lEntry,sStartGene,sStopGene,sEntryType):
    iStartCoord = 0
    iStopCoord = 0
    bStart = False
    bStop = False
    #note: last check for not bStart/bStop is in case there are multiple features with the same identifier, with different start/stop coordinates.
    #in this case only the first one should be taken.
    #why this matters: signal peptide with locus tag same as the gene on the rev. complement led to an improper determination of the stop coordinate.
    #this might still cause issues, I think, since we can't assume that the longest item is first in the genbank file.
    #I could re-sort the list, based on start position and length...mmhhh....
    for item in lEntry:
        if (item.sLocus ==sStartGene or item.sGeneName==sStartGene or item.sProduct==sStartGene or item.sProtein==sStartGene or item.sOldLocus==sStartGene) and item.sType in sEntryType and not bStart:
            iStartCoord = item.iStart
            bStart = True
        if (item.sLocus ==sStopGene or item.sGeneName==sStopGene or item.sProduct==sStopGene or item.sProtein==sStopGene or item.sOldLocus==sStopGene) and item.sType in sEntryType and not bStop:
            iStopCoord = item.iStop
            bStop = True
        '''if item.sLocus=="CD630_26050":
            print ("found it")
            print (item.sLocus,sStopGene)
            print (item.sType,bStop)  '''          
    if (not bStart) or (not bStop):
        print ("problem finding one of the genes",sStartGene,sStopGene)
        print ("exiting")
        exit(0)
    else:
        return (iStartCoord,iStopCoord)

def fill_dict(sFile):
    '''reads a random "csv" file, and parses the input
    into a dictionary. Used for reading the custom gene label and colour files
    '''
    
    inputFile = open(sFile)
    curDic = dict()
    for lines in inputFile:
        lines = lines.strip()
        if not lines:continue        
        if lines.startswith("#"):continue
        #yeah, below, that's the best I can come up with, not the smartest
        if "," in lines:sSep=","
        if ";" in lines:sSep=";"
        if "\t" in lines:sSep="\t"        
        lData = lines.split(sSep,1)
        curDic.setdefault(lData[0],lData[1])
    return curDic

def do_reverse(lEntry,iStartCoord,iStopCoord):
    for item in lEntry:
        item.iStart = iStopCoord-item.iStart+iStartCoord
        item.iStop = iStopCoord-item.iStop+iStartCoord
        if item.lIntrons:
            for i,intron in enumerate(item.lIntrons):
                #the logic for below is different than the logic for above, since here we don't care where start and stop is
                #the coordinates just need to be adjusted for the reverse system, and don't actually need to be reversed
                #there is a check for not plotting introns with negative length, so we want to have this incongruency
                item.lIntrons[i]= [iStopCoord-intron[1]+iStartCoord,iStopCoord-intron[0]+iStartCoord]
    return lEntry,iStartCoord,iStopCoord

def sanitize_output_name(sOut,sIn,sExt,sStartGene,sStopGene):
    if not sExt.startswith("."):
        sExt = "."+sExt
    if not sOut:
        sOut = sIn+"."+sStartGene+"_"+sStopGene
    return sOut,sExt

def read_input_file(sFile):
    ''' reads the file with the genbank files,
    start and stop genes and if it is reverse or not
    '''
    
    lIn = []
    lStartGene = []
    lStopGene = []
    lRev = []
    inputFile = open(sFile)
    for lines in inputFile:
        lines = lines.strip()
        if not lines:continue
        if lines.startswith("#"):continue
        lines = lines.strip("\r")
        if "," in lines:sSep=","
        if ";" in lines:sSep=";"
        if "\t" in lines:sSep="\t"
        lData  = lines.split(sSep)
        if len(lData)!=4:
            print (lines+" in "+sFile+" is missing an entry, please fix")
            exit(0)
        sFile = lData[0].strip()
        if not os.path.isfile(sFile):
            print ("could not find "+sFile+", exiting")
            exit(0)
        lIn.append(sFile)
        lStartGene.append(lData[1].strip())
        lStopGene.append(lData[2].strip())
        lRev.append(lData[3].strip())
    return lIn,lStartGene,lStopGene,lRev

def write_args(args,sFile):
    '''writes the args to a log file.
    Just in case the user wants to re-produce a plot with similar parameters.
    '''
    
    sOut = sFile+".last_input_parameters.txt"
    outputFile = open(sOut,"w")
    dicItem = vars(args)
    for key,value in dicItem.items():
        outputFile.write(key+"\t"+str(value)+"\n")
    outputFile.close()

def assign_parameters(args):
    #I could make these 4 lists into a class, but mmhh...
    #does not really feel necessary right now
    lIn = []
    lStartGene = []
    lStopGene = []
    lRev = []    
    if args.input:
        for lists in args.input:
            if not os.path.isfile(lists[0]):
                print ("could not find "+lists[0]+", exiting")
                exit(0)
            lIn.append(lists[0])
            lStartGene.append(lists[1])
            lStopGene.append(lists[2])
            lRev.append(lists[3])
    elif args.input_file:
        lIn,lStartGene,lStopGene,lRev = read_input_file(args.input_file)

    sEntryType = args.entry_type
    sLabel = args.label
    sLabelPos = args.label_location
    sOut = args.output
    sColorFile = args.color_file
    sNameFile = args.name_file
    iScale = args.scale
    iRotation = args.label_rotation
    sExt = args.file_extension
    iDistOffset = args.label_offset
    iSizeText = args.font_size
    bCoord = args.deactivate_coordinates
    fThick = args.arrow_thickness
    sOut,sExt = sanitize_output_name(sOut,lIn[0],sExt,lStartGene[0],lStopGene[0])

    return lIn,lStartGene,lStopGene,lRev,sEntryType,sLabel,sLabelPos,sOut,sColorFile,sNameFile,iScale,iRotation,sOut,sExt,iDistOffset,iSizeText,bCoord,fThick

def do_processing(args):
    '''starts the main processing
    Assigns the input variables, reads the input, sets the scale of the plot
    calls the plotting and saves the file    
    '''

    print ("version: "+VERSION)
    lIn,lStartGene,lStopGene,lRev,sEntryType,sLabel,sLabelPos,sOut,sColorFile,sNameFile,iScale,iRotation,sOut,sExt,iDistOffset,iSizeText,bCoord,fThick = assign_parameters(args)
    write_args(args,sOut)

    dicColor = dict()
    dicNames = dict()
    if sColorFile:
        dicColor = fill_dict(sColorFile)
    if sNameFile:
        dicNames = fill_dict(sNameFile)

    ###it's kinda inefficient to read all the genbanks here and further below
    ###but the logic is simply easier like this
    if not iScale:
        for i,sIn in enumerate(lIn):
            sStartGene = lStartGene[i]
            sStopGene = lStopGene[i]
            lEntry,sOrgName = read_genbank(sIn,dicColor)
            iStartCoord,iStopCoord = get_start_stop_coords(lEntry,sStartGene,sStopGene,sEntryType)            
            iDif = iStopCoord-iStartCoord
            if iDif >iScale:
                iScale = iDif
            #print (sIn,iStartCoord,iStopCoord,iDif)
    print ("longest stretch of DNA is: ",iScale)
    fFactor = 0.00320 #empirical, tested... well, guessed
    xLen = fFactor*iScale
    fig = plt.figure(figsize=(xLen,(xLen/6)*len(lIn))) #as above, just guessed around

    for i,sIn in enumerate(lIn):
        sStartGene = lStartGene[i]
        sStopGene = lStopGene[i]
        lEntry,sOrgName = read_genbank(sIn,dicColor)
        iStartCoord,iStopCoord = get_start_stop_coords(lEntry,sStartGene,sStopGene,sEntryType)
        if lRev[i]=="reverse":
            lEntry,iStartCoord,iStopCoord = do_reverse(lEntry,iStartCoord,iStopCoord)
        plot = fig.add_subplot(len(lIn),1,1+i)
        ax= plt.gca()
        make_plot(lEntry,lRev[i],iScale,sLabel,sLabelPos,iRotation,sEntryType, sStartGene,sStopGene,sOut,sExt,iStartCoord,iStopCoord,iDistOffset,iSizeText,dicNames,sOrgName,ax,bCoord,fThick)
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig(sOut+sExt)
    print (sOut+" plotted succesfully")
    plt.cla()
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", action='append', nargs=4,metavar=('genbank_file','start_gene','stop_gene',"reverse"),help="Genbank file and start and stop gene and if it should be plotted forward or reverse. Can be used multiple times to plot multiple gbk files and/or gene ranges below each other. Used -input 1.gbk locus_tag_1 locus_tag_6 forward -input 2.gbk locus_tag_8 locus_tag_12 reverse")
    parser.add_argument("--input_file",type=str,help="csv file with location of genbank files, start and stop genes and reverse/forward instructions")
    parser.add_argument('--entry_type',  nargs='+', default=["CDS","rRNA","tRNA"],help="what should be printed? CDS? genes? Any valid genbank entries will do. Default: CDS rRNA tRNA")
    parser.add_argument("--output", help="If none is given, the output will be written to inputfile+startgene+stopgene.png. Please give the complete path otherwise, this script is not smart",
                    type=str,nargs='?')
    parser.add_argument('--file_extension', help="file extension for the produced plot, default png",choices=['svg', 'png', 'pdf','jpg','ps','eps'],
                        default="png",nargs='?')
    parser.add_argument("--label_rotation", help="degree in which the labels over the genes should be rotated. default 0",
                    type=int,default=0,nargs='?')
    parser.add_argument("--font_size", help="font size used in the plots. Default 18",
                    type=int,default=18,nargs='?')
    parser.add_argument("--scale", help="the size in bp of the picture. This can be useful if multiple pictures are made, and sizes of the arrows and fonts need to be compatible. Set it to e.g. 10000 if all your plots are range 6000-9000 bp. For a single plot, all gene ranges will be scaled according to the longest, unless it is overridden like this",
                    type=int,default=0,nargs='?')
    parser.add_argument("--color_file", help="A csv file can be given with locus or gene or gene product to hex color. Tab, comma and semicolon are valid separators. Partial lists can be given. Otherwise default color coding will take place",
                    type=str,nargs='?')
    parser.add_argument("--name_file", help="A csv file can be given with locus or gene or gene product to a custom name. Tab, comma and semicolon are valid separators. Partial lists can be given. Otherwise the specified naming scheme will be used",
                    type=str,nargs='?')
    parser.add_argument("--label_offset", help="Specify if you want to have the labels higher or lower. Default 0. The height of the picture ranges from 0.1 to -0.1",
                    type=float,default=0,nargs='?')
    parser.add_argument('--label_location', help="Should the label be above or below the genes?",choices=['Up', 'Down'],
                        default="Up",nargs='?')
    parser.add_argument('--label', help="What should be printed as label? gene_name will be in italics. This can be overriden with a csv file",choices=['gene_name','locus','product','locus+product','locus+gene_name','gene_name+product'],
                        default="gene_name",nargs='?')
    parser.add_argument("--arrow_thickness", help="Factor by which the arrow should be fattened. 1.5 means the arrow will be 50 percent thicker",
                    type=float,default=1,nargs='?')    
    parser.add_argument('--deactivate_coordinates', help="Deactivate the display of genomic coordinates to the left and right of the first and last gene",action='store_false')
    parser.add_argument('-v','--version', action='store_true')
    args = parser.parse_args()
    
    if args.version:
        print ("version: "+VERSION)
        exit(1)
    if args.input or args.input_file:
        do_processing(args)
    else:
        print ("you need to provide input file and start/stop genes either via --input or --input_file")
        parser.print_help()
        exit(0)

    



        


