---
title: Getting Started - DNAnexus Developer Workshop
keywords: getting_started
sidebar: tutorial_sidebar
permalink: index.html
---

## Workshop Goals

In this workshop, we're going to walk through how to implement tools on the
DNAnexus platform. We'll identify best practices, patterns to avoid, and
provide some templates and code snippets useful for future work.

At the end of the workshop, we'll have built an executable that runs a
somewhat complex pipeline efficiently and securely on DNAnexus.

## Prerequisites

To get the most out of this workshop, you'll need a few things:

1. **DNAnexus account** - During the workshop we'll be building and running some
tools on the DNAnexus platform.

2. **Access to the workshop project** - Example input and reference files are on
the DNAnexus platform. You can use them to run and debug your tools.

3. **Access to the workshop Github repo** - The code from various stages of our
tool building process are available on Github.

4. **Basic experience with DNAnexus** - This isn't strictly necessary, but we'll
skip over some of the basics of working on DNAnexus in order to focus on
building applets.

## Workshop Task

We are going to implement an RNA-seq analysis pipeline recently published in
[*Nature Protocols*](http://www.nature.com/nprot/journal/v11/n9/full/nprot.2016.095.html):

{% include image.html file="nature_protocols_abstract.png" url="http://www.nature.com/nprot/journal/v11/n9/full/nprot.2016.095.html" alt="Nature Protocols abstract" caption="" %}

The pipeline comprises three tools:
[HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) for read alignment,
[StringTie](http://www.ccb.jhu.edu/software/stringtie/index.shtml) for
transcript assembly and quantification, and
[Ballgown](https://github.com/alyssafrazee/ballgown) for differential
expression analysis. The authors have provided a schematic describing how the
different tools fit together within their recommended pipeline:

{% include image.html file="DE_Pipeline.png" caption="Analysis schematic" %}

## Workshop Structure

We'll build up the executable that runs the pipeline iteratively. As we expand
on each component, we'll address issues that arise when developing on DNAnexus
and find solutions to common

Throughout the

{% include links.html %}
