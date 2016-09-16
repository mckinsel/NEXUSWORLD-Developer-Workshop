---
title: Organizing a Pipeline
keywords: getting_started
sidebar: tutorial_sidebar
permalink: tutorial_organizing_a_pipeline.html
---

## Breaking a Pipeline into Applets

Since we'll want our RNA-seq pipeline to be reproducible and easy to run, we're
going to build an applet. But how many applets? And if there's more than one
applet, how should work be divided among them?

Recall the schematic of the pipeline we're implementing:

{% include image.html file="DE_Pipeline.png" caption="Analysis schematic" %}

There are five different executions of a command line tool, and some of those
executions are run per sample. On one extreme, we could run everything in one
applet, which would mean it would all run in a single remote computer. On the
other extreme, we could break everything into separate applets, so that each
sample/tool combination gets its own applet. There are trade-offs with each
approach.

### Applet Trade-offs

Better recovery and restartability
: When work is split up among different applets, it's easier to restart the
  pipeline from the middle if one step should fail. If everything is
  run within one applet, the pipeline will have to begin from the start.

Better tuning of resources
: Some steps of a pipeline may require a lot of CPU or a lot of memory while
  others do not. Running steps in separate applets allows the selection of
  instance types most appropriate for each step.

Better modularity
: This can be double-edged sword, but more applets means that it's easier to
  reuse and test components of your pipeline. Taken too far, however, it can
  make it laborious to accomplish simple tasks.

Worse overhead
: Each applet runs on a separate worker, and each worker takes some time to
  get started. This can range from around 15 seconds to several tens of minutes.
  Breaking up a pipeline into many small, quick-running applets can mean that
  most of the resources dedicated to running the pipeline are actually used
  starting up workers.

Worse file I/O
: Each applet runs on a separate worker, and workers don't have a shared
  filesystem between them. So, input and reference files generally have to be
  downloaded to each worker. Having many, small applets can increase running
  time by increasing the time spent transferring files to and from workers.

{% include antipattern.html content="Sometimes users coming from environments like SGE intuit that very finely breaking down tasks into many applets will improve execution speed on DNAnexus. But unlike may grid engine environments, the overhead to starting a job on DNAnexus is rather high. So, this will actually greatly slow down execution time. If job running times in a long-running, multi-step pipeline are less than ten minutes or so, consider coalescing steps into fewer applets." %}
## RNA-seq Pipeline Organization

For our RNA-seq pipeline, we'll break it up into three applets: one that runs
HISAT2, one that runs the three StringTie steps, and one that runs Ballgown.

{% include image.html file="DE_Pipeline_broken_up.png" caption="Division into applets" %}

Ballgown is an R application, so separating it lets us address its resource
needs separately. Splitting up HISAT2 and StringTie gives us a nice checkpoint
between to steps that are likely compute-intensive. And combining the three
StringTie steps together means that we don't have to transfer assembled
transcripts and merged transcripts back and forth from platform to worker.
