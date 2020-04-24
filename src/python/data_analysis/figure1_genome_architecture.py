import glob
from PIL import Image, ImageDraw, ImageFont
import yaml

#Setting defaults
gene1_x_left = 150
gene1_x_right = 650
gene1_y_bottom = 275
gene1_y_top = 225

gene2_x_left = 750
gene2_x_right = 1250
gene2_y_bottom = 275
gene2_y_top = 225

gene3_x_left = 1350
gene3_x_right = 1850
gene3_y_bottom = 275
gene3_y_top = 225
genome_end = 1949

element_y_top = 205
element_y_bottom = 295
y_mid_area = 250
offset = 50

def create_template_genome():
    """
    Genome architecture of the template genome is produced.
    """

    im = Image.new('RGB', (2000, 500), "WHITE")
    draw = ImageDraw.Draw(im)

    #Draw genes

    #Gene 1 with blue color
    draw.rectangle((gene1_x_left, gene1_y_bottom, gene1_x_right, gene1_y_top), fill=(86, 180, 233), outline=(86, 180, 233))
    #Gene 2 with orange color
    draw.rectangle((gene2_x_left, gene2_y_bottom, gene2_x_right, gene2_y_top), fill=(230, 159, 0), outline=(230, 159, 0))
    #Gene 3 with purple color
    draw.rectangle((gene3_x_left, gene3_y_bottom, gene3_x_right, gene3_y_top), fill=(204, 121, 167), outline=(204, 121, 167))

    #Draw promoters (arrows)
    draw.line(((100, element_y_bottom), (100, element_y_top)), fill=(0, 0, 0), width=6)
    draw.line(((97, element_y_top), (130, element_y_top)), fill=(0, 0, 0), width=6)
    draw.polygon(((150, 205), (130, 195), (130, 215)), fill=(0, 0, 0))

    #Draw intergenic regions (line segments)
    draw.line(((103, y_mid_area), (gene1_x_left-1, y_mid_area)), fill=(220, 220, 220), width=2)
    draw.line(((gene1_x_right+1, y_mid_area), (gene2_x_left-1, y_mid_area)), fill=(220, 220, 220), width=2)
    draw.line(((gene2_x_right+1, y_mid_area), (gene3_x_left-1, y_mid_area)), fill=(220, 220, 220), width=2)
    draw.line(((gene3_x_right+1, y_mid_area), (genome_end, y_mid_area)), fill=(220, 220, 220), width=2)

    #Insert text
    font = ImageFont.truetype("./fonts/arial-italic.ttf", 20)
    draw.text((((gene1_x_right+gene1_x_left)/2) - 25, gene1_y_bottom + 5), 'gene 1', fill=(0, 0, 0), font=font)
    draw.text((((gene2_x_right+gene2_x_left)/2) - 25, gene2_y_bottom + 5), 'gene 2', fill=(0, 0, 0), font=font)
    draw.text((((gene3_x_right+gene3_x_left)/2) - 25, gene3_y_bottom + 5), 'gene 3', fill=(0, 0, 0), font=font)

    #Save image
    im.save('../../../data/figures/figure1/figure1_template_genome_arch.png', quality=100)

def create_genome_architecture(gene_file):
    """
    Create genome architecture based on input file.
    """

    #Setting defaults
    gene1_x_left = 150
    gene1_x_right = 650
    gene1_y_bottom = 275
    gene1_y_top = 225

    gene2_x_left = 750
    gene2_x_right = 1250
    gene2_y_bottom = 275
    gene2_y_top = 225

    gene3_x_left = 1350
    gene3_x_right = 1850
    gene3_y_bottom = 275
    gene3_y_top = 225
    genome_end = 1949

    element_y_top = 205
    element_y_bottom = 295
    y_mid_area = 250
    offset = 50

    im = Image.new('RGB', (3000, 500), "WHITE")
    draw = ImageDraw.Draw(im)

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
                elif '2' in promoter:
                    region_two[promoter] = gene_file[promoter]['start']
        if gene != 0:
            if gene_file[terminator]['start'] > 0:
                if '1' in terminator:
                    region_one[terminator] = gene_file[terminator]['start']
                elif '2' in terminator:
                    region_two[terminator] = gene_file[terminator]['start']
                elif '3' in terminator:
                    region_three[terminator] = gene_file[terminator]['start']
        if gene != gene_file['num_genes']:
            if gene_file[rnase]['start'] > 0:
                if '0' in rnase:
                    region_zero[rnase] = gene_file[rnase]['start']
                elif '1' in rnase:
                    region_one[rnase] = gene_file[rnase]['start']
                elif '2' in rnase:
                    rregion_two[rnase] = gene_file[rnase]['start']
    #Draw elements depending on their location
    if region_zero != {}:
        #Draw rnase
        draw.line(((155, element_y_top), (gene1_x_left+49, gene1_y_bottom+20)), fill=(0, 0, 0), width=6)
        draw.line(((gene1_x_left+49, element_y_top), (155, gene1_y_bottom+20)), fill=(0, 0, 0), width=6)
        #Offset all the elements to the right of it to make space for drawing
        gene1_x_left += offset
        gene1_x_right += offset
        gene2_x_left += offset
        gene2_x_right += offset
        gene3_x_left += offset
        gene3_x_right += offset
        genome_end += offset
    if region_one != {}:
        pass
    if region_two != {}:
        pass
    if region_three != {}:
        draw.line(((gene3_x_right+103, element_y_bottom), (gene3_x_right+103, element_y_top)), fill=(0, 0, 0), width=6)
        draw.line(((gene3_x_right+70, element_y_top), (gene3_x_right+133, element_y_top)), fill=(0, 0, 0), width=6)

    #Draw genes

    #Gene 1 with blue color
    draw.rectangle((gene1_x_left, gene1_y_bottom, gene1_x_right, gene1_y_top), fill=(86, 180, 233), outline=(86, 180, 233))
    #Gene 2 with orange color
    draw.rectangle((gene2_x_left, gene2_y_bottom, gene2_x_right, gene2_y_top), fill=(230, 159, 0), outline=(230, 159, 0))
    #Gene 3 with purple color
    draw.rectangle((gene3_x_left, gene3_y_bottom, gene3_x_right, gene3_y_top), fill=(204, 121, 167), outline=(204, 121, 167))

    #Draw promoters (arrows)
    draw.line(((100, element_y_bottom), (100, element_y_top)), fill=(0, 0, 0), width=6)
    draw.line(((97, element_y_top), (130, element_y_top)), fill=(0, 0, 0), width=6)
    draw.polygon(((150, 205), (130, 195), (130, 215)), fill=(0, 0, 0))

    #Draw intergenic regions (line segments)
    draw.line(((103, y_mid_area), (gene1_x_left-1, y_mid_area)), fill=(220, 220, 220), width=2)
    draw.line(((gene1_x_right+1, y_mid_area), (gene2_x_left-1, y_mid_area)), fill=(220, 220, 220), width=2)
    draw.line(((gene2_x_right+1, y_mid_area), (gene3_x_left-1, y_mid_area)), fill=(220, 220, 220), width=2)
    draw.line(((gene3_x_right+1, y_mid_area), (genome_end, y_mid_area)), fill=(220, 220, 220), width=2)

    #Insert text
    font = ImageFont.truetype("./fonts/arial-italic.ttf", 20)
    draw.text((((gene1_x_right+gene1_x_left)/2) - 25, gene1_y_bottom + 5), 'gene 1', fill=(0, 0, 0), font=font)
    draw.text((((gene2_x_right+gene2_x_left)/2) - 25, gene2_y_bottom + 5), 'gene 2', fill=(0, 0, 0), font=font)
    draw.text((((gene3_x_right+gene3_x_left)/2) - 25, gene3_y_bottom + 5), 'gene 3', fill=(0, 0, 0), font=font)

    #Save image
    im.save('../../../data/figures/figure1/figure1_intermed_genome_arch.png', quality=100)


def main():

    create_template_genome()
    with open('../../../results/2020_4_16/positive_control_nf0.7_rep0_nmut5/gene_0.yml', 'r') as gene_parameters:
        gene_file = yaml.safe_load(gene_parameters)
    create_genome_architecture(gene_file)



main()
