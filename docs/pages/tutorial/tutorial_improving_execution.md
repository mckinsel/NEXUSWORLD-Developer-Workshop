---
title: Improving Execution
keywords: hisat2
sidebar: tutorial_sidebar
permalink: tutorial_improving_execution.html
---

## Instance Type

### Choosing an instance type

## Ubuntu 14.04

The default OS used by DNAnexus workers is Ubuntu 12.04. This was released in April 2012 and
will soon reach EOL.

```json
{
  "runSpec": {
    "distribution": "Ubuntu",
    "release": "14.04"
    ...
  }
}
```

## Timeout Policy


```json
{
  "runSpec": {
    "timeoutPolicy": {"*": {"hours": 12}}
    ...
  }
}
```
