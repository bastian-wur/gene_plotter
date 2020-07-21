# gene_plotter
Gene plotter is a tool to plot/draw ranges of genes from one or multiple genbank files.
It allows custom labels, colouring, and a range of other options, mostly dervied from the underlying matplotlib library.

An example can be seen below:<br/>

![example output of gene_plotter](/example/Clostridium_difficile_paloc.png)

This plot was generated with the following command:<br/>
python3 gene_plotter.py --input_file example/plot.csv --color_file example/colour_file.txt --file_extension png --name_file example/name_file.txt --font_size 36 --label_rotation 25 --label_offset -0.01 --arrow_thickness 2<br/>

This tool has been tested with a range of bacterial genomes and a bit of yeast, so probably not all bugs have been found yet. Feedback is welcome.

## Requirements
This tool is written in Python3.<br/>
It requires matplotlib.<br/>

## How to use
Minimum use:<br/>
python3 gene_plotter.py --input_file /home/plot1.csv<br/>
OR <br/>
python3 gene_plotter.py --input genbank1.gb gene_start gene_end forward -input genbank2.gb gene_start2 gene_end2 reverse<br/>

Longer example:<br/>
python3 gene_plotter.py --input_file plot1.csv --color_file /home/colour_file.txt --file_extension svg --name_file /home/name_file.txt --font_size 42

## Options
* --input genbank_file start_gene stop_gene reverse<br/>
                        Genbank file and start and stop gene and if it should
                        be plotted forward or reverse. Can be used multiple
                        times to plot multiple gbk files and/or gene ranges
                        below each other. Used -i 1.gbk locus_tag_1
                        locus_tag_6 forward -i 2.gbk locus_tag_8 locus_tag_12
                        reverse.
                        As start/stop identifiers the following can be used:
                        locus_tag, old_locus_tag, gene, product, protein_id
*  --input_file INPUT_FILE<br/>
                        csv file with location of genbank files, start and
                        stop genes and reverse/forward instructions. Lines starting with # will be ignored.
*  --entry_type ENTRY_TYPE [ENTRY_TYPE ...]<br/>
                        what should be printed? CDS? genes? Any valid genbank
                        entries will do. Default: CDS rRNA tRNA
*  --output [OUTPUT]<br/>
                        If none is given, the output will be written to
                        inputfile1+startgene+stopgene.png. Please give the
                        complete path otherwise, this script is not smart
*  --file_extension [{svg,png,pdf,jpg,ps,eps}]<br/>
                        file extension for the produced plot, default png
*  --label_rotation [LABEL_ROTATION]<br/>
                        degree in which the labels over the genes should be
                        rotated. default 0
*  --font_size [FONT_SIZE]<br/>
                        font size used in the plots. Default 18
*  --scale [SCALE]<br/>
                        the size in bp of the picture. This can be useful if
                        multiple pictures are made, and sizes of the arrows
                        and fonts need to be compatible. Set it to e.g. 10000
                        if all your plots are range 6000-9000 bp. For a single
                        plot, all gene ranges will be scaled according to the
                        longest, unless it is overridden like this. Above 50000 
                        this might become not so very useful
*  --color_file [COLOR_FILE]<br/>
                        A csv file can be given with locus or gene or gene
                        product to hex color. Tab, comma and semicolon are
                        valid separators. Partial lists can be given.
                        Otherwise default color coding will take place.
                        Matplotlib defaul colors will be accepted (black, white, red, green, etc, see https://matplotlib.org/3.1.0/gallery/color/named_colors.html),
                        otherwise hex colours, #ADD8E6 or #FFFFFF.  Lines starting with # will be ignored.
*  --name_file [NAME_FILE]<br/>
                        A csv file can be given with locus or gene or gene
                        product to a custom name. Tab, comma and semicolon are
                        valid separators. Partial lists can be given.
                        Otherwise the specified naming scheme will be used. 
                        To make items italics, you will need to use
                        e.g. $\it{tetR}$. Lines starting with # will be ignored.
*  --label_offset [LABEL_OFFSET]<br/>
                        Specify if you want to have the labels higher or
                        lower. Default 0
*  --label_location [{Up,Down}]<br/>
                        Should the label be above or below the genes?
*  --label [{gene_name,locus,product,locus+product,locus+gene_name,gene_name+product}]<br/>
                        What should be printed as label? gene_name will be in
                        italics. This can be overriden with a csv file (no
                        styling though). Default is gene_name, otherwise locus_tag
*  --deactivate_coordinates<br/>
                        Deactivate the display of genomic coordinates to the
                        left and right of the first and last gene
*  --arrow_thickness<br/>
                        Factor by which the arrow should be fattened. 1.5 means the arrow will be 50% thicker                      
*  -v, --version<br/>                        

## What does it not do
This tool ONLY draws the genes into a plot.<br/>
It does not show in any form homology between genes, you need to know which genes belong to each other.<br/>
While this could be interesting to add, it would add considerable overhead, so the chance of implementation is low.<br/>

## Current bugs
1. italics for organism names: There is a heuristic to determine what should and should not be in italics, and it should catch most cases, but it is not perfect. To fix wrong italics in Inkscape manually, select the character, then go object -> transform -> skew -> horizontal, and use as value 11
1. complicated genbank files: to avoid dependencies, the parser for genbank files is simple and self made. It might break with more comlicated genbank entries.
1. overlapping features with the same identifer: If there are features with the same identifier, which overlap (e.g. locus tag used for both CDS and for signal peptide), then this will lead to weird results.
1. There is some weird scaling issue with the height of displayed introns. This does not become apparent for most cases, but sometimes there is a bit of an offset.
<br/>
Please still report bugs, in case you encounter them from genbank files downloaded from a widespread source, or if some obvious check for strain identifiers is missing. 

## Features not implemented (yet)
1. vertical alignment: You might want to align the different plot parts to e.g. a gene in the middle. This is currently not implemented, since this probably complicates the logic of the plotting considerably. This feature will probably not be implemented.
1. proper error handling: Some errors in the input files will not be caught right now

## Citation
Please cite this github repository if you use gene_plotter to produce any figures in your paper.<br/>
A publication in maybe the journal of open source software JOSS is planned for the future.<br/>

## License
This software is distributed under the GPLv3.<br/>

