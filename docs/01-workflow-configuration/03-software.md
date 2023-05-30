---
title: Software
permalink: workflow-configuration/software
layout: page
nav_order: 2
parent: Workflow configuration
---

# Software
{: .no_toc }

A breakdown of available software-environment management options and a list of pinned software versions.
{: .fs-6 .fw-300 }

Software
{: .label }



<details markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>

---
## Docker Environments
A docker image has been created using [Alpine Linux](), a bare-bones version of linux optimized for overall install size and ideal OS for lightweight Docker images. This image is available on [DockerHub here]().
<!-- TODO: Provide Link -->
<!-- TODO: Provide Link -->

## Conda Environments
A theoretical conda environment exists in the project, however it is not actively tested to due available work capacity. You are welcome to use it at your own risk. As capacity becomes availabel this will be revisited.

## Virtual Environments
Pinned package version files for `pip` are available be default (maintained via the `pip-tools` packages `pip-compile --generate-hashes` command). These files can be used to reliably recreate a consistent software environment with pre-existing tools as follows:

**Linux:**
```bash
# IF no venv folder present THEN
# CREATE virtual environment
python3 -m venv venv

# ACTIVATE virtual environment 'venv'
source ./venv/bin/activate

# (venv) INSTALL dependancies
pip install -r requirements.deb.txt
```

{: .normal-title}
> Virtual environments
>
> Users are encouraged to use Python virtual environments to maintain a clean, project-specific environment. To create one, you can run the `python -m venv <folder-to-create>` command to create a vritual environment. To activate this environment, use the `source <path-to<folder-to-create>>/bin/activate`. Now, any installation commands will install to this blank environment.

{: .normal-title }
> Windows compatability
>
> Virtual environments are supported on windows, however, the `source` command is linux-specific. Users can use an IDE which is capable of auto-detecting and auto-enabling virtual environments (such as [Visual Studio Code]()) for ease-of-use.