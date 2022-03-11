---
title: "Biochar Stability Dashboard"
description: "
<p>
<b>Part 1 / 3   -   Exploring the data</b>
<p/>
<p>
This first dashboard allows to explore the biochar incubation data and the associated metadata, which were collected from the litterature and compiled in a single consistent format. The result of many hours of work.
</p>

"
bg_image: "/img/lin-field.JPG"
draft: false

type : dashboard
layout : dashbord

menu:
  main:
    parent: "More"
    name: "Biochar stability dashboard"
    weight: 2

# Page-specific JavaScript & CSS, from static folder #ESA
js : ["bokehapp/bokeh-2.4.3.min", "bokehapp/bokeh-widgets-2.4.3.min", "bokehapp/bokeh-tables-2.4.3.min", "bokehapp/bokeh-gl-2.4.3.min", "bokehapp/bokeh-mathjax-2.4.3.min", ]
css : ["bokehapp/bokehapp-custom"]

---
<div class="col-md-12">
{{% notice warning %}}
Dashboard under construction 
{{% /notice %}}
</div>

<!-- MAIN DASHBOARD -->
<section class="">
  <div class="container-fluid">
 
<div class="row">
  <div class="col-md-6 col-md-push-1">

## **Biochar incubation data and the associated metadata**
  
The two graphs below present the actual biochar incubation data, for all observations in the data set, both as cumulative biochar carbon loss (decay, in %) and as daily decay rates (% biochar C loss per day). The linked table provides some information about the experimental conditions (what we call metadata).

*Figures are interactive: hover on the lines, click on the lines, select rows from the table, to highlight specific curves.*


  </div>
</div>

<div class="row">
<div id="component_1" style="height:90vh; margin:0 5% 1% 5%; max-width: 90%;" >
    <br />
    <br />
</div>
</div>

</div>
</section>

<!-- DATA DOWNLOAD  -->
<section class="section-sm section-dark">
<div class="container">
<div class="row">
  <div class="col-md-9 col-md-push-1">

## **Download all the data**

All the data shown above (and the code to deal with it) is available for download in various format:
* Get it from *GitHub* (data + code)
* Get it in *.xlsx* format (data only)
* Get it in *.json* format (data only)

## **Scroll down for more**

Below, more interactive figures and text analyses are provided. They reproduce the results from our pre-print article published here *insert link*.

  </div>
</div>
</div>
</section>

<!-- FURTHER FIGURES -->
<section class="">
<div class="container-fluid">

<div class="row">
  <div class="col-md-6 col-md-push-1">

## **What biomass types were used?**

  </div>
</div>

<div class="row">
<div class="col-md-7">
<div id="component_2" style="max-height:90vh; margin:0 5% 1% 5%;" >
    <br />
</div>
</div>

<div class="col-md-3">
<div id="component_text_2">

### **Most biochars were derived from wood**

Most observations were from wood-derived biochars (74 observations). Other biomass types used were crop residues (29, mainly corncob, rice straw or wheat straw), grass (18, mainly miscanthus), manure (4) and biosolids (1). 



### **Some feedstocks are preferred for practical reasons**

Biomass with a C4 carbon fixation pathway (e.g., maize, miscanthus, sugarcane) has been commonly used as it is a cheaper and easier way to perform isotopic analyses than with labelled biomass.



### **Wood-derived biochars have much lower ash content than other feedstocks**

A key difference between the biochar from different feedstocks is the ash content. 

</div>
</div>
</div>

</div>
</section>

<!-- FURTHER FIGURES -->
<section class="">
<div class="container-fluid">

<div class="row">
  <div class="col-md-6 col-md-push-1">

## **What pyrolysis conditions types were used?**

  </div>
</div>

<div class="row">
<div class="col-md-7">
<div id="component_3" style="max-height:90vh; margin:0 5% 1% 5%;" >
    <br />
</div>
</div>

<div class="col-md-3">
<div id="component_text_3">

### **Most biochars were from slow pyrolysis**

The vast majority of observations used biochar produced via slow pyrolysis (119 observations), while only a handful of observations are available for gasification (4), hydrothermal carbonisation (2), or flash pyrolysis (1).


### **Most biochars were produced at lab-scale**

In most cases, the pyrolysis was performed in laboratory bench scale reactors, although some commercial pyrolysis reactors have also been used (mainly when using a biomass feedstock with a C4 carbon fixating pathway, for practical reasons).

### **Pyrolysis temperature, heating rate, and residence time**

Pyrolysis conditions are usually described using the highest treatment temperature (HHT, °C), the heating rate (HR, °C/min) and the residence time (RT, min) at the HHT. The HHT is nearly always reported in incubation studies, although it is recognised that HHT is difficult to determine in large scale (commercial) reactors due to non homogenous heat distribution. Likewise, HR and RT are mostly reported for laboratory produced biochar where these parameters are well controlled. A common possible source of misinterpretation, for the RT, is that some studies report an overall RT inside the reactor (from start to end of pyrolysis process, possibly including drying time) while other report the RT at HHT. Here, only RT at HHT were included in the database leading to some missing data.

The range of HHT investigated spans from, with 200 to 800°C, with most values in the range of 400-600°C. One experiment used a gasification biochar produced at a much higher temperature of 1200°C. HR spanned between 2.5°C/min and up to 51°C/min. Finally, RT is mostly short, between 40 minutes and 3 hours, although the effect of longer residence time was tested in some cases, up to 72 hours.

From this exploration, it can be drawn that there is a gap in coverage of biochar produced at higher temperatures (> 700°C) and that coverage of the lower range of temperature could be improved as it is likely to be associated with lower stability and non-linearity. Likewise, conditions other than slow pyrolysis should also better covered (at least, in the literature performing incubation with isotopic analysis), including hydro-thermal carbonisation, gasification, and fast/flash pyrolysis.

</div>
</div>
</div>

</div>
</section>

<section>
</section>

<div id="bokeh_component_scripts">
</div>