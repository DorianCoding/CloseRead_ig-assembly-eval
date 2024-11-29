import re
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import warnings
import pandas as pd
from fpdf import FPDF
import cairosvg
from IPython.display import SVG, display
import math
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from fpdf import FPDF
from PIL import Image
from intervaltree import Interval, IntervalTree

warnings.filterwarnings('ignore')

class PDF(FPDF):
    def header(self):
        self.set_font('Arial', 'B', 12)

    def rotate(self, angle, x=None, y=None):
        if x is None:
            x = self.x
        if y is None:
            y = self.y
        self._out('Q')  # Restore previous graphics state
        self._out('q')  # Save graphics state
        if angle != 0:
            angle *= -1  # Invert angle to match expected rotation direction
            c = math.cos(math.radians(angle))
            s = math.sin(math.radians(angle))
            cx = x * self.k
            cy = (self.h - y) * self.k
            self._out(f'1 0 0 1 {cx:.2f} {cy:.2f} cm')  # Translate
            self._out(f'{c:.5f} {s:.5f} {-s:.5f} {c:.5f} 0 0 cm')  # Rotate
            self._out(f'1 0 0 1 {-cx:.2f} {-cy:.2f} cm')  # Translate back

    def stop_rotation(self):
        self._out('Q')

def make_pdf(both,image_files, output_filename, species, CommonName, LatinName, hapkind, Source, Link, overlapx=0.87, overlapy=0.01, scale_top=0.75, scale_bottom=0.3):

    # Open the first image to get dimensions
    first_image = Image.open(image_files[0])
    width, height = first_image.size
    first_width = width * 0.8
    first_height = height* 0.8
    pagewidth = width*1.15 if both else width*1.6
    pageheight = height*2.1 if both else height*1.6
    # Create instance of FPDF class, adjust the page size based on expected content
    pdf = PDF(unit="pt", format=[pagewidth, pageheight])  # Increased page height to fit vertical stacking
    w = pdf.w - 2 * pdf.l_margin
    pdf.add_page()
    pdf.set_font("Arial", 'I', 20)  # Bold Arial, 24 pt
    pdf.image("https://www.imgt.org/images/logoIMGT_medium.png",0,None,0,70,"PNG","https://imgt.org")
    pdf.set_fill_color(0)
    pdf.cell(0,20,"Results from CloseRead software (GNUv3 license)",0,0,"L",False,"https://github.com/IMGT-CNRS/CloseRead_ig-assembly-eval")
    pdf.image("https://github.githubassets.com/assets/GitHub-Mark-ea2971cee799.png",w - pdf.r_margin,None,0,60,"PNG","https://github.com/IMGT-CNRS/CloseRead_ig-assembly-eval")
    pdf.set_font("Arial", 'B', 32)  # Bold Arial, 24 pt
    pdf.set_text_color(0, 0, 0)  # Black text
    # Add the first image at full scale
    size = 850 if both else 1050
    pdf.image(image_files[0], size, 50, first_width, first_height)

    # Open the second image to get dimensions and calculate new dimensions based on the scale
    second_image = Image.open(image_files[1])
    second_width, second_height = second_image.size
    new_width_top = second_width * scale_top
    new_height_top = second_height * scale_top

    # Add the second image, resized and positioned with overlap
    overlap_x = width * overlapx
    overlap_y = height * overlapy
    pdf.image(image_files[1], x=overlap_x+size, y=overlap_y+50, w=new_width_top, h=new_height_top)

    initial_y_position_below = first_height + 100
    # Add images below x and y, vertically stacked
    current_y = initial_y_position_below
    img = Image.open(image_files[2])
    orig_width, orig_height = img.size
    scaled_width = width*0.933 if both else width*1.2 # Scale each image to fit the page width
    scaled_height = orig_height * (scaled_width / orig_width)  # Maintain aspect ratio
    pdf.image(image_files[2], x=800, y=current_y, w=scaled_width, h=scaled_height)
    current_y += scaled_height  # Update the y position for the next image
    max=6 if both else 4
    for image_path in image_files[3:max]:
        if image_path == image_files[3] or image_path == image_files[5]:
            img = Image.open(image_path)
            orig_width, orig_height = img.size
            scaled_width = width * 0.85 if both else width * 1.1  # Scale each image to fit the page width
            scaled_height = orig_height * (scaled_width / orig_width)  # Maintain aspect ratio
        else:
            img = Image.open(image_path)
            orig_width, orig_height = img.size
            scaled_width = width * 0.925 if both else width * 1.1  # Scale each image to fit the page width
            scaled_height = orig_height * (scaled_width / orig_width)  # Maintain aspect ratio
        if both:
            pdf.image(image_path, x=800, y=current_y, w=scaled_width, h=scaled_height)
        else:
            pdf.image(image_path, x=775, y=current_y+50, w=scaled_width, h=scaled_height)
        current_y += scaled_height  # Update the y position for the next image
    pdf.set_font("Arial", '', 36)  
    title = f"Species ID: {species}"
    commonName = f"Common Name: {CommonName}"
    sciName = f"Scientific Name: {LatinName}"
    hap = f"Assembly Type: {hapkind}"
    source = f"Data Source: {Source}"
    title_width = pdf.get_string_width(title) + 6
    pdf.set_y(200)
    pdf.cell(125)
    pdf.cell(title_width, 20, title, 0, 1, 'L')
    pdf.set_y(300)
    pdf.cell(125)
    pdf.cell(title_width, 20, commonName, 0, 1, 'L')
    pdf.set_y(400)
    pdf.cell(125)
    pdf.cell(title_width, 20, sciName, 0, 1, 'L')
    pdf.set_y(500)
    pdf.cell(125)
    pdf.cell(title_width, 20, hap, 0, 1, 'L')
    pdf.set_y(600)
    pdf.cell(125)
    pdf.cell(title_width, 20, source, 0, 1, "L", False,Link)
    # Save the result to a PDF file
    pdf.set_subject("CloseRead results")
    pdf.set_creator("CloseRead software")
    pdf.output(output_filename)
