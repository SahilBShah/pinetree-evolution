import glob
from PIL import Image, ImageDraw, ImageFont
import yaml

#Setting defaults
gene1_x_left = 108
gene1_x_right = 508
gene1_y_bottom = 105
gene1_y_top = 55

gene2_x_left = 608
gene2_x_right = 1008
gene2_y_bottom = 105
gene2_y_top = 55

gene3_x_left = 1108
gene3_x_right = 1508
gene3_y_bottom = 105
gene3_y_top = 55

genome_end = 1607

element_y_top = 14
y_mid_area = 80
offset = 50

def create_template_genome():
    """
    Genome architecture of the template genome is produced.
    """

    im = Image.new('RGB', (1610, 200), "WHITE")
    draw = ImageDraw.Draw(im)

    #Draw genes

    #Gene 1 with blue color
    draw.rectangle((gene1_x_left, gene1_y_bottom, gene1_x_right, gene1_y_top), fill=(86, 180, 233), outline=(86, 180, 233))
    #Gene 2 with orange color
    draw.rectangle((gene2_x_left, gene2_y_bottom, gene2_x_right, gene2_y_top), fill=(230, 159, 0), outline=(230, 159, 0))
    #Gene 3 with purple color
    draw.rectangle((gene3_x_left, gene3_y_bottom, gene3_x_right, gene3_y_top), fill=(204, 121, 167), outline=(204, 121, 167))

    #Draw promoters (arrows)
    draw.line(((5, y_mid_area+12), (5, element_y_top)), fill=(0, 0, 0), width=6)
    draw.line(((2, element_y_top), (38, element_y_top)), fill=(0, 0, 0), width=6)
    draw.polygon(((68, element_y_top), (38, element_y_top-15), (38, element_y_top+15)), fill=(0, 0, 0))

    #Draw intergenic regions (line segments)
    draw.line(((8, y_mid_area), (gene1_x_left-1, y_mid_area)), fill=(192, 192, 192), width=25)
    draw.line(((gene1_x_right+1, y_mid_area), (gene2_x_left-1, y_mid_area)), fill=(192, 192, 192), width=25)
    draw.line(((gene2_x_right+1, y_mid_area), (gene3_x_left-1, y_mid_area)), fill=(192, 192, 192), width=25)
    draw.line(((gene3_x_right+1, y_mid_area), (genome_end, y_mid_area)), fill=(192, 192, 192), width=25)

    #Insert text
    font = ImageFont.truetype("./fonts/arial-italic.ttf", 80)
    draw.text((((gene1_x_right+gene1_x_left)/2) - 130, gene1_y_bottom + 5), 'gene X', fill=(0, 0, 0), font=font)
    draw.text((((gene2_x_right+gene2_x_left)/2) - 130, gene2_y_bottom + 5), 'gene Y', fill=(0, 0, 0), font=font)
    draw.text((((gene3_x_right+gene3_x_left)/2) - 125, gene3_y_bottom + 5), 'gene Z', fill=(0, 0, 0), font=font)

    #Save image
    im.save('../../../figure_output/figure1/figure1_template_genome_arch.png', quality=100)

def create_ten_gene_genome():

    #Setting defaults
    gene1_x_left = 108
    gene1_x_right = 508
    gene1_y_bottom = 105
    gene1_y_top = 55

    gene2_x_left = 608
    gene2_x_right = 1008
    gene2_y_bottom = 105
    gene2_y_top = 55

    gene3_x_left = 1108
    gene3_x_right = 1508
    gene3_y_bottom = 105
    gene3_y_top = 55

    gene4_x_left = 1608
    gene4_x_right = 2008
    gene4_y_bottom = 105
    gene4_y_top = 55

    gene5_x_left = 2108
    gene5_x_right = 2508
    gene5_y_bottom = 105
    gene5_y_top = 55

    gene6_x_left = 2608
    gene6_x_right = 3008
    gene6_y_bottom = 105
    gene6_y_top = 55

    gene7_x_left = 3108
    gene7_x_right = 3508
    gene7_y_bottom = 105
    gene7_y_top = 55

    gene8_x_left = 3608
    gene8_x_right = 4008
    gene8_y_bottom = 105
    gene8_y_top = 55

    gene9_x_left = 4108
    gene9_x_right = 4508
    gene9_y_bottom = 105
    gene9_y_top = 55

    gene10_x_left = 4608
    gene10_x_right = 5008
    gene10_y_bottom = 105
    gene10_y_top = 55

    genome_end = 5107

    element_y_top = 14
    y_mid_area = 80

    for count in range(1, 3):

        im = Image.new('RGB', (5110, 200), "WHITE")
        draw = ImageDraw.Draw(im)

        #Draw genes

        #Gene 1 with blue color
        draw.rectangle((gene1_x_left, gene1_y_bottom, gene1_x_right, gene1_y_top), fill=(0, 70, 139), outline=(0, 70, 139))
        #Gene 2 with orange color
        draw.rectangle((gene2_x_left, gene2_y_bottom, gene2_x_right, gene2_y_top), fill=(0, 153, 180), outline=(0, 153, 180))
        #Gene 3 with purple color
        draw.rectangle((gene3_x_left, gene3_y_bottom, gene3_x_right, gene3_y_top), fill=(146, 94, 159), outline=(146, 94, 159))
        #Gene 4 with a color
        draw.rectangle((gene4_x_left, gene4_y_bottom, gene4_x_right, gene4_y_top), fill=(253, 175, 145), outline=(253, 175, 145))
        #Gene 5 with b color
        draw.rectangle((gene5_x_left, gene5_y_bottom, gene5_x_right, gene5_y_top), fill=(173, 0, 42), outline=(173, 0, 42))
        #Gene 6 with c color
        draw.rectangle((gene6_x_left, gene6_y_bottom, gene6_x_right, gene6_y_top), fill=(173, 182, 182), outline=(173, 182, 182))
        #Gene 7 with d color
        draw.rectangle((gene7_x_left, gene7_y_bottom, gene7_x_right, gene7_y_top), fill=(27, 25, 25), outline=(27, 25, 25))
        #Gene 8 with e color
        draw.rectangle((gene8_x_left, gene8_y_bottom, gene8_x_right, gene8_y_top), fill=(237, 0, 0), outline=(237, 0, 0))
        #Gene 9 with f color
        draw.rectangle((gene9_x_left, gene9_y_bottom, gene9_x_right, gene9_y_top), fill=(2, 75, 48), outline=(2, 75, 48))
        #Gene 10 with g color
        draw.rectangle((gene10_x_left, gene10_y_bottom, gene10_x_right, gene10_y_top), fill=(66, 181, 64), outline=(66, 181, 64))

        #Draw promoters (arrows)
        draw.line(((5, y_mid_area+12), (5, element_y_top)), fill=(0, 0, 0), width=6)
        draw.line(((2, element_y_top), (38, element_y_top)), fill=(0, 0, 0), width=6)
        draw.polygon(((68, element_y_top), (38, element_y_top-15), (38, element_y_top+15)), fill=(0, 0, 0))

        if count == 2:

            #Draw terminators
            #Terminator 2
            draw.line(((gene2_x_right+33, y_mid_area), (gene2_x_right+33, element_y_top)), fill=(0, 0, 0), width=6)
            draw.line(((gene2_x_right, element_y_top), (gene2_x_right+63, element_y_top)), fill=(0, 0, 0), width=6)
            #Terminator 3
            draw.line(((gene3_x_right+33, y_mid_area), (gene3_x_right+33, element_y_top)), fill=(0, 0, 0), width=6)
            draw.line(((gene3_x_right, element_y_top), (gene3_x_right+63, element_y_top)), fill=(0, 0, 0), width=6)
            #Terminator 4
            draw.line(((gene4_x_right+33, y_mid_area), (gene4_x_right+33, element_y_top)), fill=(0, 0, 0), width=6)
            draw.line(((gene4_x_right, element_y_top), (gene4_x_right+63, element_y_top)), fill=(0, 0, 0), width=6)
            #Terminator 6
            draw.line(((gene6_x_right+33, y_mid_area), (gene6_x_right+33, element_y_top)), fill=(0, 0, 0), width=6)
            draw.line(((gene6_x_right, element_y_top), (gene6_x_right+63, element_y_top)), fill=(0, 0, 0), width=6)
            #Terminator 7
            draw.line(((gene7_x_right+33, y_mid_area), (gene7_x_right+33, element_y_top)), fill=(0, 0, 0), width=6)
            draw.line(((gene7_x_right, element_y_top), (gene7_x_right+63, element_y_top)), fill=(0, 0, 0), width=6)
            #Terminator 8
            draw.line(((gene8_x_right+33, y_mid_area), (gene8_x_right+33, element_y_top)), fill=(0, 0, 0), width=6)
            draw.line(((gene8_x_right, element_y_top), (gene8_x_right+63, element_y_top)), fill=(0, 0, 0), width=6)

            #Draw Promoters
            #Promoter 6
            draw.line(((gene6_x_right+85, y_mid_area), (gene6_x_right+85, element_y_top)), fill=(0, 0, 0), width=6)
            draw.line(((gene6_x_right+82, element_y_top), (gene6_x_right+118, element_y_top)), fill=(0, 0, 0), width=6)
            draw.polygon(((gene6_x_right+148, element_y_top), (gene6_x_right+118, element_y_top-15), (gene6_x_right+118, element_y_top+15)), fill=(0, 0, 0))

        #Draw intergenic regions (line segments)
        draw.line(((8, y_mid_area), (gene1_x_left-1, y_mid_area)), fill=(192,192,192), width=25)
        draw.line(((gene1_x_right+1, y_mid_area), (gene2_x_left-1, y_mid_area)), fill=(192,192,192), width=25)
        draw.line(((gene2_x_right+1, y_mid_area), (gene3_x_left-1, y_mid_area)), fill=(192,192,192), width=25)
        draw.line(((gene3_x_right+1, y_mid_area), (gene4_x_left-1, y_mid_area)), fill=(192,192,192), width=25)
        draw.line(((gene4_x_right+1, y_mid_area), (gene5_x_left-1, y_mid_area)), fill=(192,192,192), width=25)
        draw.line(((gene5_x_right+1, y_mid_area), (gene6_x_left-1, y_mid_area)), fill=(192,192,192), width=25)
        draw.line(((gene6_x_right+1, y_mid_area), (gene7_x_left-1, y_mid_area)), fill=(192,192,192), width=25)
        draw.line(((gene7_x_right+1, y_mid_area), (gene8_x_left-1, y_mid_area)), fill=(192,192,192), width=25)
        draw.line(((gene8_x_right+1, y_mid_area), (gene9_x_left-1, y_mid_area)), fill=(192,192,192), width=25)
        draw.line(((gene9_x_right+1, y_mid_area), (gene10_x_left-1, y_mid_area)), fill=(192,192,192), width=25)
        draw.line(((gene10_x_right+1, y_mid_area), (genome_end, y_mid_area)), fill=(192,192,192), width=25)

        #Insert text
        font = ImageFont.truetype("./fonts/arial-italic.ttf", 80)
        draw.text((((gene1_x_right+gene1_x_left)/2) - 130, gene1_y_bottom + 5), 'gene A', fill=(0, 0, 0), font=font)
        draw.text((((gene2_x_right+gene2_x_left)/2) - 130, gene2_y_bottom + 5), 'gene B', fill=(0, 0, 0), font=font)
        draw.text((((gene3_x_right+gene3_x_left)/2) - 125, gene3_y_bottom + 5), 'gene C', fill=(0, 0, 0), font=font)
        draw.text((((gene4_x_right+gene4_x_left)/2) - 125, gene4_y_bottom + 5), 'gene D', fill=(0, 0, 0), font=font)
        draw.text((((gene5_x_right+gene5_x_left)/2) - 125, gene5_y_bottom + 5), 'gene E', fill=(0, 0, 0), font=font)
        draw.text((((gene6_x_right+gene6_x_left)/2) - 125, gene6_y_bottom + 5), 'gene F', fill=(0, 0, 0), font=font)
        draw.text((((gene7_x_right+gene7_x_left)/2) - 125, gene7_y_bottom + 5), 'gene G', fill=(0, 0, 0), font=font)
        draw.text((((gene8_x_right+gene8_x_left)/2) - 125, gene8_y_bottom + 5), 'gene H', fill=(0, 0, 0), font=font)
        draw.text((((gene9_x_right+gene9_x_left)/2) - 125, gene9_y_bottom + 5), 'gene I', fill=(0, 0, 0), font=font)
        draw.text((((gene10_x_right+gene10_x_left)/2) - 125, gene10_y_bottom + 5), 'gene J', fill=(0, 0, 0), font=font)

        if count == 2:
            #Draw RNAses
            #RNAse 9
            draw.line(((gene9_x_right+10, gene9_y_top-7), (gene9_x_right+58, gene9_y_bottom+7)), fill=(0, 0, 0), width=6)
            draw.line(((gene9_x_right+58, gene9_y_top-7), (gene9_x_right+10, gene9_y_bottom+7)), fill=(0, 0, 0), width=6)

        #Save image
        if count == 2:
            im.save('../../../figure_output/figure5/figure5_genome_arch_final.png', quality=100)
        else:
            im.save('../../../figure_output/figure5/figure5_genome_arch_starting.png', quality=100)

def create_genome_architecture(gene_file, output_dir):
    """
    Create genome architecture based on input file.
    """

    #Setting defaults
    gene1_x_left = 108
    gene1_x_right = 508
    gene1_y_bottom = 105
    gene1_y_top = 55

    gene2_x_left = 608
    gene2_x_right = 1008
    gene2_y_bottom = 105
    gene2_y_top = 55

    gene3_x_left = 1108
    gene3_x_right = 1508
    gene3_y_bottom = 105
    gene3_y_top = 55

    genome_end = 1607

    element_y_top = 14
    y_mid_area = 80

    prom_offset = 10
    term_offset = 35
    rnase_offset = 20
    total_offset = 1610

    region_zero = {}
    region_one = {}
    region_two = {}
    region_three = {}

    #Enumerate elements on genome
    for gene in range(gene_file['num_genes']+1):
        promoter = 'promoter_{}'.format(gene)
        terminator = 'terminator_{}'.format(gene)
        rnase = 'rnase_{}'.format(gene)
        if gene != gene_file['num_genes'] and gene != 0:
            if gene_file[promoter]['start'] > 0:
                if '1' in promoter:
                    region_one[promoter] = gene_file[promoter]['start']
                    total_offset+=40
                elif '2' in promoter:
                    region_two[promoter] = gene_file[promoter]['start']
                    total_offset+=40
        if gene != 0:
            if gene_file[terminator]['start'] > 0:
                if '1' in terminator:
                    region_one[terminator] = gene_file[terminator]['start']
                    total_offset+=40
                elif '2' in terminator:
                    region_two[terminator] = gene_file[terminator]['start']
                    total_offset+=40
                elif '3' in terminator:
                    region_three[terminator] = gene_file[terminator]['start']
                    total_offset+=40
        if gene != gene_file['num_genes']:
            if gene_file[rnase]['start'] > 0:
                if '0' in rnase:
                    region_zero[rnase] = gene_file[rnase]['start']
                    total_offset+=40
                elif '1' in rnase:
                    region_one[rnase] = gene_file[rnase]['start']
                    total_offset+=40
                elif '2' in rnase:
                    region_two[rnase] = gene_file[rnase]['start']
                    total_offset+=40

    im = Image.new('RGB', (total_offset, 200), "WHITE")
    draw = ImageDraw.Draw(im)

    #Variables needed to determine the number of elements to draw
    region_one_total = len(region_one)
    region_two_total = len(region_two)
    increment_one = 0
    increment_two = 0
    rnase_flags = [False, False, False, False, False, False, False, False]


    #Draw elements depending on their location
    if region_zero != {}:
        #Draw rnase
        draw.line(((65, gene1_y_top-7), (gene1_x_left+5, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
        draw.line(((gene1_x_left+5, gene1_y_top-7), (65, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
        #Offset all the elements to the right of it to make space for drawing
        gene1_x_left += 15
        gene1_x_right += 15
        gene2_x_left += 15
        gene2_x_right += 15
        gene3_x_left += 15
        gene3_x_right += 15
        genome_end += 15
        rnase_flags[0] = True
    if region_one != {}:
        while region_one != {}:
            if increment_one == 0:
                min_element = min(region_one.keys(), key=(lambda k: region_one[k]))
                del region_one[min_element]
                if 'promoter' in min_element:
                    draw.line(((gene1_x_right+5, y_mid_area), (gene1_x_right+5, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.line(((gene1_x_right+2, element_y_top), (gene1_x_right+38, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.polygon(((gene1_x_right+68, element_y_top), (gene1_x_right+38, element_y_top-15), (gene1_x_right+38, element_y_top+15)), fill=(0, 0, 0))
                    gene2_x_left += prom_offset
                    gene2_x_right += prom_offset
                    gene3_x_left += prom_offset
                    gene3_x_right += prom_offset
                    genome_end += prom_offset
                elif 'terminator' in min_element:
                    draw.line(((gene1_x_right+33, y_mid_area), (gene1_x_right+33, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.line(((gene1_x_right, element_y_top), (gene1_x_right+63, element_y_top)), fill=(0, 0, 0), width=6)
                    gene2_x_left += term_offset
                    gene2_x_right += term_offset
                    gene3_x_left += term_offset
                    gene3_x_right += term_offset
                    genome_end += term_offset
                elif 'rnase' in min_element:
                    draw.line(((gene1_x_right+10, gene1_y_top-7), (gene1_x_right+58, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
                    draw.line(((gene1_x_right+58, gene1_y_top-7), (gene1_x_right+10, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
                    gene2_x_left += rnase_offset
                    gene2_x_right += rnase_offset
                    gene3_x_left += rnase_offset
                    gene3_x_right += rnase_offset
                    genome_end += rnase_offset
                    rnase_flags[1] = True
            if increment_one == 1:
                min_element = min(region_one.keys(), key=(lambda k: region_one[k]))
                del region_one[min_element]
                if 'promoter' in min_element:
                    draw.line(((gene1_x_right+85, y_mid_area), (gene1_x_right+85, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.line(((gene1_x_right+82, element_y_top), (gene1_x_right+118, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.polygon(((gene1_x_right+148, element_y_top), (gene1_x_right+118, element_y_top-15), (gene1_x_right+118, element_y_top+15)), fill=(0, 0, 0))
                    gene2_x_left += prom_offset
                    gene2_x_right += prom_offset
                    gene3_x_left += prom_offset
                    gene3_x_right += prom_offset
                    genome_end += prom_offset
                elif 'terminator' in min_element:
                    draw.line(((gene1_x_right+107, y_mid_area), (gene1_x_right+107, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.line(((gene1_x_right+74, element_y_top), (gene1_x_right+139, element_y_top)), fill=(0, 0, 0), width=6)
                    gene2_x_left += term_offset
                    gene2_x_right += term_offset
                    gene3_x_left += term_offset
                    gene3_x_right += term_offset
                    genome_end += term_offset
                elif 'rnase' in min_element:
                    draw.line(((gene1_x_right+80, gene1_y_top-7), (gene1_x_right+128, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
                    draw.line(((gene1_x_right+128, gene1_y_top-7), (gene1_x_right+80, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
                    gene2_x_left += rnase_offset
                    gene2_x_right += rnase_offset
                    gene3_x_left += rnase_offset
                    gene3_x_right += rnase_offset
                    genome_end += rnase_offset+5
                    rnase_flags[2] = True
            if increment_one == 2:
                min_element = min(region_one.keys(), key=(lambda k: region_one[k]))
                del region_one[min_element]
                if 'promoter' in min_element:
                    draw.line(((gene1_x_right+155, y_mid_area), (gene1_x_right+155, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.line(((gene1_x_right+152, element_y_top), (gene1_x_right+188, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.polygon(((gene1_x_right+218, element_y_top), (gene1_x_right+188, element_y_top-15), (gene1_x_right+188, element_y_top+15)), fill=(0, 0, 0))
                    gene2_x_left += prom_offset
                    gene2_x_right += prom_offset
                    gene3_x_left += prom_offset
                    gene3_x_right += prom_offset
                    genome_end += prom_offset+30
                elif 'terminator' in min_element:
                    draw.line(((gene1_x_right+196, y_mid_area), (gene1_x_right+196, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.line(((gene1_x_right+161, element_y_top), (gene1_x_right+227, element_y_top)), fill=(0, 0, 0), width=6)
                    gene2_x_left += term_offset+70
                    gene2_x_right += term_offset+70
                    gene3_x_left += term_offset+70
                    gene3_x_right += term_offset+70
                    genome_end += term_offset+40
                elif 'rnase' in min_element:
                    draw.line(((gene1_x_right+150, gene1_y_top-7), (gene1_x_right+198, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
                    draw.line(((gene1_x_right+198, gene1_y_top-7), (gene1_x_right+150, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
                    gene2_x_left += rnase_offset+40
                    gene2_x_right += rnase_offset+40
                    gene3_x_left += rnase_offset+40
                    gene3_x_right += rnase_offset+40
                    genome_end += rnase_offset+40
                    rnase_flags[3] = True
            increment_one+=1
    if region_two != {}:
        while region_two != {}:
            if increment_two == 0:
                min_element = min(region_two.keys(), key=(lambda k: region_two[k]))
                del region_two[min_element]
                if 'promoter' in min_element:
                    draw.line(((gene2_x_right+5, y_mid_area), (gene2_x_right+5, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.line(((gene2_x_right+2, element_y_top), (gene2_x_right+38, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.polygon(((gene2_x_right+68, element_y_top), (gene2_x_right+38, element_y_top-15), (gene2_x_right+38, element_y_top+15)), fill=(0, 0, 0))
                    gene3_x_left += prom_offset+5
                    gene3_x_right += prom_offset+5
                    genome_end += prom_offset+5
                elif 'terminator' in min_element:
                    draw.line(((gene2_x_right+33, y_mid_area), (gene2_x_right+33, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.line(((gene2_x_right, element_y_top), (gene2_x_right+63, element_y_top)), fill=(0, 0, 0), width=6)
                    gene3_x_left += term_offset
                    gene3_x_right += term_offset
                    genome_end += term_offset
                elif 'rnase' in min_element:
                    draw.line(((gene2_x_right+10, gene1_y_top-7), (gene2_x_right+58, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
                    draw.line(((gene2_x_right+58, gene1_y_top-7), (gene2_x_right+10, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
                    gene3_x_left += rnase_offset
                    gene3_x_right += rnase_offset
                    genome_end += rnase_offset
                    rnase_flags[4] = True
            if increment_two == 1:
                min_element = min(region_two.keys(), key=(lambda k: region_two[k]))
                del region_two[min_element]
                if 'promoter' in min_element:
                    draw.line(((gene2_x_right+85, y_mid_area), (gene2_x_right+85, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.line(((gene2_x_right+82, element_y_top), (gene2_x_right+118, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.polygon(((gene2_x_right+148, element_y_top), (gene2_x_right+118, element_y_top-15), (gene2_x_right+118, element_y_top+15)), fill=(0, 0, 0))
                    gene3_x_left += prom_offset+5
                    gene3_x_right += prom_offset+5
                    genome_end += prom_offset+5
                elif 'terminator' in min_element:
                    draw.line(((gene2_x_right+107, y_mid_area), (gene2_x_right+107, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.line(((gene2_x_right+74, element_y_top), (gene2_x_right+139, element_y_top)), fill=(0, 0, 0), width=6)
                    gene3_x_left += term_offset
                    gene3_x_right += term_offset
                    genome_end += term_offset
                elif 'rnase' in min_element:
                    draw.line(((gene2_x_right+80, gene1_y_top-7), (gene2_x_right+128, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
                    draw.line(((gene2_x_right+128, gene1_y_top-7), (gene2_x_right+80, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
                    gene3_x_left += rnase_offset
                    gene3_x_right += rnase_offset
                    genome_end += rnase_offset
                    rnase_flags[5] = True
            if increment_two == 2:
                min_element = min(region_two.keys(), key=(lambda k: region_two[k]))
                del region_two[min_element]
                if 'promoter' in min_element:
                    draw.line(((gene2_x_right+155, y_mid_area), (gene2_x_right+155, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.line(((gene2_x_right+152, element_y_top), (gene2_x_right+188, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.polygon(((gene2_x_right+218, element_y_top), (gene2_x_right+188, element_y_top-15), (gene2_x_right+188, element_y_top+15)), fill=(0, 0, 0))
                    gene3_x_left += prom_offset+60
                    gene3_x_right += prom_offset+60
                    genome_end += prom_offset+60
                elif 'terminator' in min_element:
                    draw.line(((gene2_x_right+195, y_mid_area), (gene2_x_right+195, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.line(((gene2_x_right+160, element_y_top), (gene2_x_right+226, element_y_top)), fill=(0, 0, 0), width=6)
                    gene3_x_left += term_offset+30
                    gene3_x_right += term_offset+30
                    genome_end += term_offset+30
                elif 'rnase' in min_element:
                    draw.line(((gene2_x_right+150, gene1_y_top-7), (gene2_x_right+198, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
                    draw.line(((gene2_x_right+198, gene1_y_top-7), (gene2_x_right+150, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
                    gene3_x_left += rnase_offset
                    gene3_x_right += rnase_offset
                    genome_end += rnase_offset
                    rnase_flags[6] = True
            increment_two+=1
    if region_three != {}:
        draw.line(((genome_end-2, y_mid_area), (genome_end-2, element_y_top)), fill=(0, 0, 0), width=6)
        draw.line(((genome_end-35, element_y_top), (genome_end+28, element_y_top)), fill=(0, 0, 0), width=6)

    #Draw genes

    #Gene 1 with blue color
    draw.rectangle((gene1_x_left, gene1_y_bottom, gene1_x_right, gene1_y_top), fill=(86, 180, 233), outline=(86, 180, 233))
    #Gene 2 with orange color
    draw.rectangle((gene2_x_left, gene2_y_bottom, gene2_x_right, gene2_y_top), fill=(230, 159, 0), outline=(230, 159, 0))
    #Gene 3 with purple color
    draw.rectangle((gene3_x_left, gene3_y_bottom, gene3_x_right, gene3_y_top), fill=(204, 121, 167), outline=(204, 121, 167))

    #Draw promoters (arrows)
    draw.line(((5, y_mid_area+12), (5, element_y_top)), fill=(0, 0, 0), width=6)
    draw.line(((2, element_y_top), (38, element_y_top)), fill=(0, 0, 0), width=6)
    draw.polygon(((68, element_y_top), (38, element_y_top-15), (38, element_y_top+15)), fill=(0, 0, 0))

    #Draw intergenic regions (line segments)
    draw.line(((8, y_mid_area), (gene1_x_left-1, y_mid_area)), fill=(192,192,192), width=25)
    draw.line(((gene1_x_right+1, y_mid_area), (gene2_x_left-1, y_mid_area)), fill=(192,192,192), width=25)
    draw.line(((gene2_x_right+1, y_mid_area), (gene3_x_left-1, y_mid_area)), fill=(192,192,192), width=25)
    draw.line(((gene3_x_right+1, y_mid_area), (genome_end, y_mid_area)), fill=(192,192,192), width=25)

    #Insert text
    font = ImageFont.truetype("./fonts/arial-italic.ttf", 80)
    draw.text((((gene1_x_right+gene1_x_left)/2) - 130, gene1_y_bottom + 5), 'gene X', fill=(0, 0, 0), font=font)
    draw.text((((gene2_x_right+gene2_x_left)/2) - 130, gene2_y_bottom + 5), 'gene Y', fill=(0, 0, 0), font=font)
    draw.text((((gene3_x_right+gene3_x_left)/2) - 125, gene3_y_bottom + 5), 'gene Z', fill=(0, 0, 0), font=font)

    #Redraw RNAses
    for flag in range(len(rnase_flags)):
        if flag == 0 and rnase_flags[flag] == True:
            draw.line(((65, gene1_y_top-7), (gene1_x_left-10, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
            draw.line(((gene1_x_left-10, gene1_y_top-7), (65, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
        if flag == 1 and rnase_flags[flag] == True:
            draw.line(((gene1_x_right+10, gene1_y_top-7), (gene1_x_right+58, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
            draw.line(((gene1_x_right+58, gene1_y_top-7), (gene1_x_right+10, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
        if flag == 2 and rnase_flags[flag] == True:
            draw.line(((gene1_x_right+80, gene1_y_top-7), (gene1_x_right+128, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
            draw.line(((gene1_x_right+128, gene1_y_top-7), (gene1_x_right+80, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
        if flag == 3 and rnase_flags[flag] == True:
            draw.line(((gene1_x_right+150, gene1_y_top-7), (gene1_x_right+198, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
            draw.line(((gene1_x_right+198, gene1_y_top-7), (gene1_x_right+150, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
        if flag == 4 and rnase_flags[flag] == True:
            draw.line(((gene2_x_right+10, gene1_y_top-7), (gene2_x_right+58, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
            draw.line(((gene2_x_right+58, gene1_y_top-7), (gene2_x_right+10, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
        if flag == 5 and rnase_flags[flag] == True:
            draw.line(((gene2_x_right+80, gene1_y_top-7), (gene2_x_right+128, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
            draw.line(((gene2_x_right+128, gene1_y_top-7), (gene2_x_right+80, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
        if flag == 6 and rnase_flags[flag] == True:
            draw.line(((gene2_x_right+150, gene1_y_top-7), (gene2_x_right+198, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
            draw.line(((gene2_x_right+198, gene1_y_top-7), (gene2_x_right+150, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)

    #Save image
    im.save(output_dir, quality=100)


def main():

    #Figure 1
    #Template genome architecture
    create_template_genome()
    #Target genome architecture
    with open('../../figure_data/genome_configurations/positive_control.yml', 'r') as gene_parameters:
        gene_file = yaml.safe_load(gene_parameters)
    output_dir = '../../figure_output/figure1/figure1_beg_genome_arch.png'
    create_genome_architecture(gene_file, output_dir)
    #Generation 0 genome architecture
    with open('../../manuscript_results/fig1/positive_control/rep5/gene_0.yml', 'r') as gene_parameters:
        gene_file = yaml.safe_load(gene_parameters)
    output_dir = '../../figure_output/figure1/figure1_gen0_genome_arch.png'
    create_genome_architecture(gene_file, output_dir)
    #Intermediate generation genome architecture
    with open('../../manuscript_results/fig1/positive_control/rep5/gene_488.yml', 'r') as gene_parameters:
        gene_file = yaml.safe_load(gene_parameters)
    output_dir = '../../figure_output/figure1/figure1_intermed_genome_arch.png'
    create_genome_architecture(gene_file, output_dir)
    #Final genome architecture
    with open('../../manuscript_results/fig1/positive_control/rep5/final/gene_best.yml', 'r') as gene_parameters:
        gene_file = yaml.safe_load(gene_parameters)
    output_dir = '../../figure_output/figure1/figure1_final_genome_arch.png'
    create_genome_architecture(gene_file, output_dir)

    #Figure 2
    #Final first genome architecture
    with open('../../manuscript_results/fig2/paper_data1_arrange1/rep18/final/gene_best.yml', 'r') as gene_parameters:
        gene_file = yaml.safe_load(gene_parameters)
    output_dir = '../../figure_output/figure2/figure2_final1_genome_arch.png'
    create_genome_architecture(gene_file, output_dir)
    #Final second genome architecture
    with open('../../manuscript_results/fig2/paper_data1_arrange5/rep8/final/gene_best.yml', 'r') as gene_parameters:
        gene_file = yaml.safe_load(gene_parameters)
    output_dir = '../../figure_output/figure2/figure2_final2_genome_arch.png'
    create_genome_architecture(gene_file, output_dir)
    #Final third genome architecture
    with open('../../manuscript_results/fig2/paper_data5_arrange1/rep20/final/gene_best.yml', 'r') as gene_parameters:
        gene_file = yaml.safe_load(gene_parameters)
    output_dir = '../../figure_output/figure2/figure2_final3_genome_arch.png'
    create_genome_architecture(gene_file, output_dir)
    #Final fourth genome architecture
    with open('../../manuscript_results/fig2/paper_data6_arrange1/rep21/final/gene_best.yml', 'r') as gene_parameters:
        gene_file = yaml.safe_load(gene_parameters)
    output_dir = '../../figure_output/figure2/figure2_final4_genome_arch.png'
    create_genome_architecture(gene_file, output_dir)
    #Final fifth genome architecture
    with open('../../manuscript_results/fig2/paper_data8_arrange5/rep8/final/gene_best.yml', 'r') as gene_parameters:
        gene_file = yaml.safe_load(gene_parameters)
    output_dir = '../../figure_output/figure2/figure2_final5_genome_arch.png'
    create_genome_architecture(gene_file, output_dir)

    #Promoter genome architecture - parameter sweep
    with open('../../figure_data/genome_configurations/promoter_test.yml', 'r') as gene_parameters:
        gene_file = yaml.safe_load(gene_parameters)
    output_dir = '../../figure_output/figure_element_strengths/promoter_genome_arch.png'
    create_genome_architecture(gene_file, output_dir)
    #Terminator genome architecture - parameter sweep
    with open('../../figure_data/genome_configurations/figure_element_strengths/terminator_test.yml', 'r') as gene_parameters:
        gene_file = yaml.safe_load(gene_parameters)
    output_dir = '../../figure_output/terminator_genome_arch.png'
    create_genome_architecture(gene_file, output_dir)
    #RNAse genome architecture - parameter sweep
    with open('../../figure_data/genome_configurations/rnase_test.yml', 'r') as gene_parameters:
        gene_file = yaml.safe_load(gene_parameters)
    output_dir = '../../figure_output/figure_element_strengths/rnase_genome_arch.png'
    create_genome_architecture(gene_file, output_dir)

    create_ten_gene_genome()

    # #Creates architectures for a specific directory
    # for i in (1, 2, 3, 4, 5, 7, 8, 9, 13, 14, 18, 19, 20, 21, 22, 23, 25, 26, 27, 28, 30, 31, 32, 33, 35, 36, 37, 38, 39, 40, 41, 42, 43, 45, 47, 48, 49, 50):
	   #  with open('../../manuscript_results/paper_data1_arrange1/rep{}/final/gene_clean.yml'.format(i), 'r') as gene_parameters:
	   #       gene_file = yaml.safe_load(gene_parameters)
	   #  output_dir = '../../../results/genome_arch1_arrange1_rep{}_clean.png'.format(i)
	   #  create_genome_architecture(gene_file, output_dir)

if __name__ == "__main__":
    main()
