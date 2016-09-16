---
title: Basic HISAT2 Applet
keywords: hisat2
sidebar: tutorial_sidebar
permalink: tutorial_basic_hisat2_applet.html
---

## Applet Components

An applet on DNAnexus has three parts:

[dxapp.json](https://wiki.dnanexus.com/dxapp.json)
:  The dxapp.json file defines the interface of the applet, and it tell the
   platform some details about how the applet should be run.

Applet script
:  Each applet has a script, written in bash or python2.7, that the platform
   runs within the worker when the applet it run. Often, this script just
   calls other scripts and compiled code.

Applet resources
:  Optionally, each applet can have resource files packaged with it. These will
   be placed on the worker before the applet script is run. There are a few ways
   to package applet resources, which will be discussed later.

## HISAT2 Applet Interface

The HISAT2 applet that we're building takes RNA-seq reads and a reference index
as inputs, and it outputs alignments. Based on the HISAT2 documentation, we'll
assume that the reads come in the form of two FASTQ files, the reference index
comes in the form of a gzipped tarball, and the output is a SAM file. With this
information, we can create the first part of our dxapp.json file:

```json
{
  "name": "basic_hisat2_applet_bash",
  "title": "Basic HISAT2 Applet (bash)",
  "summary": "Simply downloads input files, runs HISAT2, and upload the output SAM.",
  "description": "This is a starting point for a more sophisticated implementation of HISAT2 on DNAnexus.",
  "inputSpec": [
    {
      "name": "hisat2_index_targz",
      "label": "HISAT2 Index Tarball",
      "class": "file",
      "help": "Tar.gz'd HISAT2 index for a reference genome, produced using hisat2-build."
    },
    {
      "name": "mate1_fastq",
      "label": "Mate 1 FASTQ",
      "class": "file",
      "help": "FASTQ file containing mate 1 reads."
    },
    {
      "name": "mate2_fastq",
      "label": "Mate 2 FASTQ",
      "class": "file",
      "help": "FASTQ file containing mate 2 reads."
    }
  ],
  "outputSpec": [
    {
      "name": "aligned_sam",
      "label": "Aligned SAM",
      "class": "file",
      "help": "SAM file with alignments reported by HISAT2"
    }
  ],
```

Note that we've specified that each input and output is of the "file" class,
meaning that the platform is expecting a single file for each.

The remainder of the dxapp.json file depends on whether we're going to write the
applet script in python or bash:


<ul id="profileTabs" class="nav nav-tabs">
    <li class="active"><a href="#python" data-toggle="tab">python</a></li>
    <li><a href="#bash" data-toggle="tab">bash</a></li>
</ul>
<div class="tab-content">
  <div role="tabpanel" class="tab-pane active" id="python">
    <div class="language-json highlighter-rouge"><pre class="highlight"><code><span class="w">  </span><span class="s2">"runSpec"</span><span class="err">:</span><span class="w"> </span><span class="p">{</span><span class="w">
        </span><span class="nt">"file"</span><span class="p">:</span><span class="w"> </span><span class="s2">"src/script.py"</span><span class="p">,</span><span class="w">
        </span><span class="nt">"interpreter"</span><span class="p">:</span><span class="w"> </span><span class="s2">"python2.7"</span><span class="w">
      </span><span class="p">}</span><span class="w">
    </span><span class="err">}</span><span class="w">
    </span></code></pre>
    </div>
  </div>
  <div role="tabpanel" class="tab-pane" id="bash">
    <pre class="highlight"><code><span class="w">  </span><span class="s2">"runSpec"</span><span class="err">:</span><span class="w"> </span><span class="p">{</span><span class="w">
        </span><span class="nt">"file"</span><span class="p">:</span><span class="w"> </span><span class="s2">"src/script.sh"</span><span class="p">,</span><span class="w">
        </span><span class="nt">"interpreter"</span><span class="p">:</span><span class="w"> </span><span class="s2">"bash"</span><span class="w">
      </span><span class="p">}</span><span class="w">
    </span><span class="err">}</span><span class="w">
    </span></code></pre>
  </div>
</div>

## HISAT2 Applet Script

The command that the applet is supposed to run is pretty straightforward:

```shell
hisat2 --dta -x <hisat2-idx> -1 <m1> -2 <m2>  -S <hit>
```

Here ```<hisat2-idx>``` refers to the "basename" of the HISAT2 reference index,
```<m1>``` and ```<m2>``` refer to fastq files with the first and second mates
of the RNA-seq reads, and ```<hit>``` refers to the output SAM file. The ```--dta```
is from the description of the protocol.

That's all easy enough, but to be able to run that command we have to first
make the inputs available on the worker. Then, we have to upload the output
SAM file back to the platform and indicate that it should be associated with the
"aligned_sam" output that the platform is expecting.

<ul id="profileTabs" class="nav nav-tabs">
    <li class="active"><a href="#script_python" data-toggle="tab">python</a></li>
    <li><a href="#script_bash" data-toggle="tab">bash</a></li>
</ul>
<div class="tab-content">
<div role="tabpanel" class="tab-pane active" id="script_python">
<div class="language-python highlighter-rouge"><pre class="highlight"><code><span class="kn">import</span> <span class="nn">os</span>
<span class="kn">import</span> <span class="nn">subprocess</span>

<span class="kn">import</span> <span class="nn">dxpy</span>

<span class="nd">@dxpy.entry_point</span><span class="p">(</span><span class="s">"main"</span><span class="p">)</span>
<span class="k">def</span> <span class="nf">main</span><span class="p">(</span><span class="n">hisat2_index_targz</span><span class="p">,</span> <span class="n">mate1_fastq</span><span class="p">,</span> <span class="n">mate2_fastq</span><span class="p">):</span>

    <span class="c"># First, download all the input files to local storage</span>
    <span class="n">dxpy</span><span class="o">.</span><span class="n">download_dxfile</span><span class="p">(</span><span class="n">hisat2_index_targz</span><span class="p">,</span> <span class="s">"hisat2_index.tar.gz"</span><span class="p">)</span>
    <span class="n">dxpy</span><span class="o">.</span><span class="n">download_dxfile</span><span class="p">(</span><span class="n">mate1_fastq</span><span class="p">,</span> <span class="s">"mate1.fastq"</span><span class="p">)</span>
    <span class="n">dxpy</span><span class="o">.</span><span class="n">download_dxfile</span><span class="p">(</span><span class="n">mate2_fastq</span><span class="p">,</span> <span class="s">"mate2.fastq"</span><span class="p">)</span>


    <span class="c"># Second, extract the index tarball</span>
    <span class="n">proc</span> <span class="o">=</span> <span class="n">subprocess</span><span class="o">.</span><span class="n">Popen</span><span class="p">([</span><span class="s">"tar"</span><span class="p">,</span> <span class="s">"xf"</span><span class="p">,</span> <span class="s">"hisat2_index.tar.gz"</span><span class="p">])</span>
    <span class="n">proc</span><span class="o">.</span><span class="n">wait</span><span class="p">()</span>

    <span class="c"># You could also use the tarfile module for this</span>
    <span class="c"># import tarfile</span>
    <span class="c"># tar = tarfile.open("hisat2_index.tar.gz")</span>
    <span class="c"># tar.extractall()</span>
    <span class="c"># tar.close()</span>

    <span class="c"># Third, figure out what the basename of the hisat2 reference index is</span>
    <span class="c"># This depends on the basename following the pattern used in the indexes</span>
    <span class="c"># distributed by the authors of HISAT2, that is the index in grch37.tar.gz</span>
    <span class="c"># will extract to grch37/genome*</span>
    <span class="n">index_basename</span> <span class="o">=</span> <span class="n">os</span><span class="o">.</span><span class="n">path</span><span class="o">.</span><span class="n">join</span><span class="p">(</span><span class="n">dxpy</span><span class="o">.</span><span class="n">DXFile</span><span class="p">(</span><span class="n">hisat2_index_targz</span><span class="p">)</span><span class="o">.</span><span class="n">name</span><span class="p">[:</span><span class="o">-</span><span class="nb">len</span><span class="p">(</span><span class="s">".tar.gz"</span><span class="p">)],</span>
                                  <span class="s">"genome"</span><span class="p">)</span>

    <span class="c"># Prepare the hisat2 command and run it.</span>
    <span class="n">hisat2_cmd_template</span> <span class="o">=</span> <span class="p">(</span><span class="s">"hisat2 -x {index_basename} -1 {mate1_fastq} -2 {mate2_fastq} "</span>
                           <span class="s">"-S {hisat2_output_sam}"</span><span class="p">)</span>
    <span class="n">hisat2_cmd</span> <span class="o">=</span> <span class="n">hisat2_cmd_template</span><span class="o">.</span><span class="n">format</span><span class="p">(</span>
        <span class="n">index_basename</span><span class="o">=</span><span class="n">index_basename</span><span class="p">,</span>
        <span class="n">mate1_fastq</span><span class="o">=</span><span class="s">"mate1.fastq"</span><span class="p">,</span>
        <span class="n">mate2_fastq</span><span class="o">=</span><span class="s">"mate2.fastq"</span><span class="p">,</span>
        <span class="n">hisat2_output_sam</span><span class="o">=</span><span class="s">"hisat2_output.sam"</span><span class="p">)</span>
    <span class="n">subprocess</span><span class="o">.</span><span class="n">check_call</span><span class="p">(</span><span class="n">hisat2_cmd</span><span class="p">,</span> <span class="n">shell</span><span class="o">=</span><span class="bp">True</span><span class="p">)</span>

    <span class="c"># Upload the output SAM file.</span>
    <span class="n">uploaded_dxfile</span> <span class="o">=</span> <span class="n">dxpy</span><span class="o">.</span><span class="n">upload_local_file</span><span class="p">(</span><span class="s">"hisat2_output.sam"</span><span class="p">)</span>

    <span class="c"># Return the ID of the uploaded SAM file associated with the "aligned_sam"</span>
    <span class="c"># field in the outputSpec in dxapp.json.</span>
    <span class="k">return</span> <span class="p">{</span><span class="s">"aligned_sam"</span><span class="p">:</span> <span class="n">uploaded_dxfile</span><span class="o">.</span><span class="n">get_id</span><span class="p">()}</span>
</code></pre>
</div>
</div>
<div role="tabpanel" class="tab-pane" id="script_bash">
<div class="language-shell highlighter-rouge"><pre class="highlight"><code><span class="c">#!/bin/bash</span>

main<span class="o">()</span> <span class="o">{</span>

    dx download <span class="s2">"</span><span class="nv">$hisat2_index_targz</span><span class="s2">"</span> -o hisat2_index.tar.gz
    dx download <span class="s2">"</span><span class="nv">$mate1_fastq</span><span class="s2">"</span> -o mate1.fastq
    dx download <span class="s2">"</span><span class="nv">$mate2_fastq</span><span class="s2">"</span> -o mate2.fastq

    tar xf hisat2_index.tar.gz

    <span class="nv">index_filename</span><span class="o">=</span><span class="k">$(</span>dx describe --name <span class="s2">"</span><span class="nv">$hisat2_index_targz</span><span class="s2">"</span><span class="k">)</span>
    <span class="nv">index_basename</span><span class="o">=</span><span class="k">${</span><span class="nv">index_filename</span><span class="p">%.tar.gz</span><span class="k">}</span>/genome

    hisat2 -x <span class="nv">$index_basename</span> -1 mate1.fastq -2 mate2.fastq -S hisat2_output.sam

    <span class="nv">uploaded_id</span><span class="o">=</span><span class="k">$(</span>dx upload hisat2_output.sam<span class="k">)</span>

    dx-jobutil-add-output aligned_sam <span class="nv">$uploaded_id</span>
<span class="o">}</span>
</code></pre>
</div>
</div>
</div>

This code downloads each input file to a hard-coded file name on the worker. It
then extracts the tarball containing the HISAT2 index. Then, it runs the ```hisat2```
command. Finally, it uploads the output SAM file and indicates that that is the
file for the "aligned_sam" output.

The trickiest part is figuring out the basename of the index. There are a number
of HISAT2 indexes available from the [authors' website](ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data).
Each of those tarballs extracts to a directory based on the name of the tarball,
and the index files within the directory have the prefix "genome". So, in either
python or bash, the applet script works with the paths and strings to get the right
parameter to pass to ```hisat2```.
