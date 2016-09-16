---
title: Do You Need an Applet?
keywords: getting_started
sidebar: tutorial_sidebar
permalink: tutorial_do_you_need_an_applet.html
---

Before delving into the details of how to build applets and workflows on
DNAnexus, it's worth considering whether you can accomplish your goal without
one. For one-off or experimental tasks, there are ways to run code and tools
without having to build a new applet. These approached give up some
reproducibility and auditability, but they're very quick and easy to use.

## Swiss Army Knife

In some cases, you want to run a command or two on some files and get the output.
For example, you might want to merge overlapping entries in a BED file using a [bedtools](http://bedtools.readthedocs.io/en/latest/)
command. You could create an applet that has a BED file input and a merged BED
file output, or you could just type the command and run it using the Swiss Army
Knife app.

The Swiss Army Knife app takes an arbitrary list of files and a command string
as inputs. It downloads all the files, runs the command, and uploads any
new files. So, we could just add our BED file as an input and type in the
appropriate command:

{% include image.html file="swiss_army_knife_example.png" caption="" alt="Swiss Army Knife" %}

{% include tip.html content="When using Swiss Army Knife, it's often useful to refer to input files using wildcards, like \*.bed, so the command can be independent of the specific input files. It also allows you to refer to multiple files at once." %}

The Swiss Army Knife app has a number of tools built-in, so you can refer to them
in commands. The list of tools and more information are available on the [apps page](https://platform.dnanexus.com/app/swiss-army-knife).

## Cloud Workstation

In other cases, you want to work interactively on a remote computer that has
fast access to files stored on the DNAnexus platform as well as the rest of the
internet. You also may want to work on a computer with specific resources, like
32 cores or 200 GiB of memory. This can be accomplished with the Cloud
Workstation app. The Cloud Workstation app simply runs a DNAnexus worker for a
specified period of time that you can SSH into.

The Cloud Workstation app is described in the [DNAnexus Wiki](https://wiki.dnanexus.com/Developer-Tutorials/Cloud-Workstations).
