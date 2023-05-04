# -*- coding: utf-8 -*-
""" setup.py; setuptools control. """

from setuptools import setup

with open("readme","rb") as f:
   long_descr= f.read().decode("utf-8")

setup(
  name= "dds_analysis",
  packages= [
     "dds_analysis",
     "dds_analysis.script",
     "dds_analysis.script.script_high"],
  entry_points = {
     "console_scripts": ['dds_analysis = dds_analysis.dds_analysis:main']} ,
  version=1.0,
  description= "Python pipeline for Differential Methylation Region (DMR) , Deifferential Expressed Genes (DEG), and SNP block analysis",
  long_description= long_descr,
  author= "Junbai Wang",
  author_email= "junbai@gmail.com"
  )
 
