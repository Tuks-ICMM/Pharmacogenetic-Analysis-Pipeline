---
title: Configuration
permalink: documentation/configuration
layout: page
nav_order: 2
has_children: false
parent: Documentation
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
<!-- ## Docker Environments
A docker image has been created using [Alpine Linux](https://www.alpinelinux.org), a bare-bones version of linux optimized for size and efficiency, making it an ideal OS for lightweight Docker images. This image is available on [DockerHub here](https://hub.docker.com/repository/docker/graemeford/pipeline-os/general).

```bash
docker pull graemeford/pipeline-os:latest
``` -->

<!-- ## Virtual Environments
Users are encouraged to use Python virtual environments to maintain a clean, project-specific environment. To create one, you can run the `python -m venv <folder-to-create>` command to create a virtual environment. IDE's such as Visual Studio Code with the appropriate python language packs enabled will be able to auto-detect this virtual environment and activate it for you. After activation, any python-related commands will install to and reference from this virtual environment.

Pinned package version files for the python package manager, `pip`, are available to prime this environment with all necessary software. This is maintained via the `pip-tools` packages `pip-compile --generate-hashes` command.


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

> `--generate-hashes` is used to record an MD5-checksum digest of the files that form part of the packages file download. This allows pip to verify the file composition on download against these records, as a security function. -->