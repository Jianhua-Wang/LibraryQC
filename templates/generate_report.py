#!/usr/bin/env python3

import sys
import re
import numpy as np
import pandas as pd
import base64
from math import modf
import datetime
pd.set_option('max_colwidth',200)

if len(sys.argv) == 1:
    sys.argv=["generate_report.py","$ngs","$lib","$pbc","$pie","$his","$distribution","$norm","$fastq1","$fastq2","$library","$vector","$pattern_str"]

ngs = sys.argv[1]
lib = sys.argv[2]
pbc = sys.argv[3]

pie = sys.argv[4]
his = sys.argv[5]
distribution = sys.argv[6]
norm = sys.argv[7]

pattern = sys.argv[12]
fastq1 = sys.argv[8]
fastq2 = sys.argv[9]
library = sys.argv[10]
vector = sys.argv[11]

today = datetime.date.today()
sample_name = re.match(r'(.*)(_1)(.*)', fastq1.split('/')[-1]).group(1)

GEN_HTML = "LibraryQC_Report.html"

def dftohtml(file):
    if type(file) == str:
        df = pd.read_table(file)
    else:
        df = file
    a = df.to_html(index=False)
    b = a.split('\\n')
    if b[3] == '      <th>Unnamed: 0</th>':
        b[3] ='       <th></th>'
    c = ''.join(b[1:-1])
    return c

def floattoint(a):
    if modf(a)[0] == 0:
        a = str(int(a))
    else:
        a = str(round(a, 6))
    return a

def pngtobase64(png):
    f = open(png,"rb")
    base64_data = base64.b64encode(f.read()).decode("utf-8")
    return base64_data

parameter = {'Parameter': ['pattern','fastq1','fastq2','library','vector'], 'Value': [pattern,fastq1,fastq2,library,vector]}
parameter_df = pd.DataFrame(data=parameter)
ngs_df = pd.read_table(ngs)
pbc_df = pd.read_table(pbc)
ngs_df['Read count'] = ngs_df['Read count'].map(lambda x: floattoint(x))
pbc_df['Value'] = pbc_df['Value'].map(lambda x: floattoint(x))
parameter_html = dftohtml(parameter_df)
lib_html = dftohtml(lib)
ngs_html = dftohtml(ngs_df)
pbc_html = dftohtml(pbc_df)

pie_html = pngtobase64(pie)
his_html = pngtobase64(his)
distribution_html = pngtobase64(distribution)
norm_html = pngtobase64(norm)

message = """<html>
    <head>
        <meta charset="utf-8"> 
        <title> LibraryQC Report</title>
        <style type="text/css">
            @media screen {
            div.summary {
            width: 18em;
            position:fixed;
            top: 3em;
            margin:1em 0 0 1em;
            }
            
            div.main {
            display:block;
            position:absolute;
            overflow:auto;
            height:auto;
            width:auto;
            top:4.5em;
            bottom:2.3em;
            left:18em;
            right:0;
            border-left: 1px solid #CCC;
            padding:0 0 0 1em;
            background-color: white;
            z-index:1;
            }
            
            div.header {
            background-color: #EEE;
            border:0;
            margin:0;
            padding: 0.5em;
            font-size: 200%;
            font-weight: bold;
            position:fixed;
            width:100%;
            top:0;
            left:0;
            z-index:2;
            }
        
            div.footer {
            background-color: #EEE;
            border:0;
            margin:0;
            padding:0.5em;
            height: 1.3em;
            overflow:hidden;
            font-size: 100%;
            font-weight: bold;
            position:fixed;
            bottom:0;
            width:100%;
            z-index:2;
            }
            
            img.indented {
            margin-left: 3em;
            }
            }
            
            @media print {
            img {
                max-width:100% !important;
                page-break-inside: avoid;
            }
            h2, h3 {
                page-break-after: avoid;
            }
            div.header {
                background-color: #FFF;
            }
            
            }
            
            body {    
            font-family: sans-serif;   
            color: #000;   
            background-color: #FFF;
            border: 0;
            margin: 0;
            padding: 0;
            }
            
            div.header {
            border:0;
            margin:0;
            padding: 0.5em;
            font-size: 200%;
            font-weight: bold;
            width:100%;
            }    
            
            #header_title {
            display:inline-block;
            float:left;
            clear:left;
            }
            #header_filename {
            display:inline-block;
            float:right;
            clear:right;
            font-size: 50%;
            margin-right:2em;
            text-align: right;
            }
        
            div.header h3 {
            font-size: 50%;
            margin-bottom: 0;
            }
            
            div.summary ul {
            padding-left:0;
            list-style-type:none;
            }
            
            div.summary ul li img {
            margin-bottom:-0.5em;
            margin-top:0.5em;
            }
                
            div.main {
            background-color: white;
            }
                
            div.module {
            padding-bottom:1.5em;
            padding-top:1.5em;
            }
                
            div.footer {
            background-color: #EEE;
            border:0;
            margin:0;
            padding: 0.5em;
            font-size: 100%;
            font-weight: bold;
            width:100%;
            }
        
        
            a {
            color: #000080;
            }
        
            a:hover {
            color: #800000;
            }
                
            h2 {
            color: #800000;
            padding-bottom: 0;
            margin-bottom: 0;
            clear:left;
            }
        
            table.dataintable {
            margin-top:15px;
            border-collapse:collapse;
            border:1px solid #aaa;
            width:50%;
            }
            table.dataintable th {
            vertical-align:baseline;
            padding:5px 15px 5px 6px;
            background-color:#3F3F3F;
            border:1px solid #3F3F3F;
            text-align:left;
            color:#fff;
            }
            table.dataintable td {
            vertical-align:text-top;
            padding:6px 15px 6px 6px;
            border:1px solid #aaa;
            }
            table.dataintable tr:nth-child(odd) {
            background-color:#F5F5F5;
            }
            table.dataintable tr:nth-child(even) {
            background-color:#fff;
            }
        
            img {
            padding-top: 0;
            margin-top: 0;
            border-top: 0;
            }
        
            
            p {
            padding-top: 0;
            margin-top: 0;
            }
        </style>
    </head>
    <body>
        <div class="header">
            <div id="header_title">
                LibraryQC Report
            </div>
            <div id="header_filename"> """+str(today)+""" <br/> """+sample_name+"""
            </div>
        </div>        
        <div class="summary">
            <h2>Summary</h2>
            <ul>
                <li>
                    <a href="#M7">Parameter</a>
                </li>
                <li>
                    <a href="#M0">LIB composition</a>
                </li>
                <li>
                    <a href="#M1">NGS composition</a>
                </li>
                <li>
                    <a href="#M2">PBC qc</a>
                </li>
                <li>
                    <a href="#M3">01pie</a>
                </li>
                <li>
                    <a href="#M4">02his</a>
                </li>
                <li>
                    <a href="#M5">03distribution</a>
                </li>
                <li>
                    <a href="#M6">04norm</a>
                </li>
            </ul>
        </div>
        <div class="main">
            <div class="module">
            <h2 id="M7">
                Parameter
            </h2>
            <table class="dataintable">"""+parameter_html+"""</table>
            </div>
            <div class="module">
            <h2 id="M0">
                LIB composition
            </h2>
            <table class="dataintable">"""+lib_html+"""</table>
            </div>
            <div class="module">
            <h2 id="M1">
                NGS composition
            </h2>
            <table class="dataintable">"""+ngs_html+"""</table>
            </div>
            <div class="module">
            <h2 id="M2">
                PBC qc
            </h2>
            <table class="dataintable">"""+pbc_html+"""</table>
            </div>
            <div class="module">
                <h2 id="M3">
                    01pie
                </h2>
                <p>
                    <img class="indented" src="data:image/png;base64,"""+pie_html+"""\" alt="01pie" width="600" height="600"/>
                </p>
            </div>
            <div class="module">
                <h2 id="M4">
                    02his
                </h2>
                <p>
                    <img class="indented" src="data:image/png;base64,"""+his_html+"""\" alt="02his" width="600" height="600"/>
                </p>
            </div>
            <div class="module">
                <h2 id="M5">
                    03distribution
                </h2>
                <p>
                    <img class="indented" src="data:image/png;base64,"""+distribution_html+"""\" alt="03distribution" width="600" height="600"/>
                </p>
            </div>
            <div class="module">
                <h2 id="M6">
                    04norm
                </h2>
                <p>
                    <img class="indented" src="data:image/png;base64,"""+norm_html+"""\" alt="04norm" width="600" height="600"/>
                </p>
            </div>
        </div>
        <div class="footer">
            Produced by Jianhua Wang
        </div>
    </body>
</html>"""

f = open(GEN_HTML,'w')
f.write(message)
f.close()