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
    draw.line(((5, y_mid_area), (5, element_y_top)), fill=(0, 0, 0), width=6)
    draw.line(((2, element_y_top), (38, element_y_top)), fill=(0, 0, 0), width=6)
    draw.polygon(((68, element_y_top), (38, element_y_top-15), (38, element_y_top+15)), fill=(0, 0, 0))

    #Draw intergenic regions (line segments)
    draw.line(((8, y_mid_area), (gene1_x_left-1, y_mid_area)), fill=(220, 220, 220), width=5)
    draw.line(((gene1_x_right+1, y_mid_area), (gene2_x_left-1, y_mid_area)), fill=(220, 220, 220), width=5)
    draw.line(((gene2_x_right+1, y_mid_area), (gene3_x_left-1, y_mid_area)), fill=(220, 220, 220), width=5)
    draw.line(((gene3_x_right+1, y_mid_area), (genome_end, y_mid_area)), fill=(220, 220, 220), width=5)

    #Insert text
    font = ImageFont.truetype("./fonts/arial-italic.ttf", 80)
    draw.text((((gene1_x_right+gene1_x_left)/2) - 130, gene1_y_bottom + 5), 'gene X', fill=(0, 0, 0), font=font)
    draw.text((((gene2_x_right+gene2_x_left)/2) - 130, gene2_y_bottom + 5), 'gene Y', fill=(0, 0, 0), font=font)
    draw.text((((gene3_x_right+gene3_x_left)/2) - 125, gene3_y_bottom + 5), 'gene Z', fill=(0, 0, 0), font=font)

    #Save image
    im.save('../../../data/figures/figure1/figure1_template_genome_arch.png', quality=100)

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
                    draw.line(((gene1_x_right+166, y_mid_area), (gene1_x_right+166, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.line(((gene1_x_right+131, element_y_top), (gene1_x_right+197, element_y_top)), fill=(0, 0, 0), width=6)
                    gene2_x_left += term_offset
                    gene2_x_right += term_offset
                    gene3_x_left += term_offset
                    gene3_x_right += term_offset
                    genome_end += term_offset+40
                elif 'rnase' in min_element:
                    draw.line(((gene1_x_right+150, gene1_y_top-7), (gene1_x_right+198, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
                    draw.line(((gene1_x_right+198, gene1_y_top-7), (gene1_x_right+150, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
                    gene2_x_left += rnase_offset
                    gene2_x_right += rnase_offset
                    gene3_x_left += rnase_offset
                    gene3_x_right += rnase_offset
                    genome_end += rnase_offset+40
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
            if increment_two == 2:
                min_element = min(region_two.keys(), key=(lambda k: region_two[k]))
                del region_two[min_element]
                if 'promoter' in min_element:
                    draw.line(((gene2_x_right+155, y_mid_area), (gene2_x_right+155, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.line(((gene2_x_right+152, element_y_top), (gene2_x_right+188, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.polygon(((gene2_x_right+218, element_y_top), (gene2_x_right+188, element_y_top-15), (gene2_x_right+188, element_y_top+15)), fill=(0, 0, 0))
                    gene3_x_left += prom_offset+5
                    gene3_x_right += prom_offset+5
                    genome_end += prom_offset+5
                elif 'terminator' in min_element:
                    draw.line(((gene2_x_right+166, y_mid_area), (gene2_x_right+166, element_y_top)), fill=(0, 0, 0), width=6)
                    draw.line(((gene2_x_right+131, element_y_top), (gene2_x_right+197, element_y_top)), fill=(0, 0, 0), width=6)
                    gene3_x_left += term_offset
                    gene3_x_right += term_offset
                    genome_end += term_offset
                elif 'rnase' in min_element:
                    draw.line(((gene2_x_right+150, gene1_y_top-7), (gene2_x_right+198, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
                    draw.line(((gene2_x_right+198, gene1_y_top-7), (gene2_x_right+150, gene1_y_bottom+7)), fill=(0, 0, 0), width=6)
                    gene3_x_left += rnase_offset
                    gene3_x_right += rnase_offset
                    genome_end += rnase_offset
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
    draw.line(((5, y_mid_area), (5, element_y_top)), fill=(0, 0, 0), width=6)
    draw.line(((2, element_y_top), (38, element_y_top)), fill=(0, 0, 0), width=6)
    draw.polygon(((68, element_y_top), (38, element_y_top-15), (38, element_y_top+15)), fill=(0, 0, 0))

    #Draw intergenic regions (line segments)
    draw.line(((8, y_mid_area), (gene1_x_left-1, y_mid_area)), fill=(220, 220, 220), width=5)
    draw.line(((gene1_x_right+1, y_mid_area), (gene2_x_left-1, y_mid_area)), fill=(220, 220, 220), width=5)
    draw.line(((gene2_x_right+1, y_mid_area), (gene3_x_left-1, y_mid_area)), fill=(220, 220, 220), width=5)
    draw.line(((gene3_x_right+1, y_mid_area), (genome_end, y_mid_area)), fill=(220, 220, 220), width=5)

    #Insert text
    font = ImageFont.truetype("./fonts/arial-italic.ttf", 80)
    draw.text((((gene1_x_right+gene1_x_left)/2) - 130, gene1_y_bottom + 5), 'gene X', fill=(0, 0, 0), font=font)
    draw.text((((gene2_x_right+gene2_x_left)/2) - 130, gene2_y_bottom + 5), 'gene Y', fill=(0, 0, 0), font=font)
    draw.text((((gene3_x_right+gene3_x_left)/2) - 125, gene3_y_bottom + 5), 'gene Z', fill=(0, 0, 0), font=font)

    #Save image
    im.save(output_dir, quality=100)


def main():

    #Template genome architecture
    create_template_genome()
    #Target genome architecture
    with open('../../../data/gene_parameters/positive_control.yml', 'r') as gene_parameters:
        gene_file = yaml.safe_load(gene_parameters)
    output_dir = '../../../data/figures/figure1/figure1_beg_genome_arch.png'
    create_genome_architecture(gene_file, output_dir)
    #Generation 0 genome architecture
    with open('../../../results/2020_5_6/positive_control_nf1.0_rep0_nmut5/gene_0.yml', 'r') as gene_parameters:
        gene_file = yaml.safe_load(gene_parameters)
    output_dir = '../../../data/figures/figure1/figure1_gen0_genome_arch.png'
    create_genome_architecture(gene_file, output_dir)
    #Intermediate generation genome architecture
    with open('../../../results/2020_5_6/positive_control_nf1.0_rep0_nmut5/gene_701.yml', 'r') as gene_parameters:
        gene_file = yaml.safe_load(gene_parameters)
    output_dir = '../../../data/figures/figure1/figure1_intermed_genome_arch.png'
    create_genome_architecture(gene_file, output_dir)
    #Final genome architecture
    with open('../../../results/2020_5_6/positive_control_nf1.0_rep0_nmut5/final/gene_best.yml', 'r') as gene_parameters:
        gene_file = yaml.safe_load(gene_parameters)
    output_dir = '../../../data/figures/figure1/figure1_final_genome_arch.png'
    create_genome_architecture(gene_file, output_dir)


main()
