---
title: Configuration
permalink: workflow/configuration
layout: page
nav_order: 1
has_children: false
parent: Workflow
---

# Configuration
{: .no_toc }

How to configure the workflow for an analysis
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

## `config` folder

To perform an analysis with this workflow, users will need to configure the workflow. This involves providing environment-related information like output locations, as well as analysis settings like reference population selection. This information is all declared and stored using the `config/manifest.json` file.

<h3>The <code>manifest.json</code> file</h3>

This file is responsible for declaring all information relating to the analysis and serves as the central point of contact between the workflow runtime and your input data. It is also used to configure and synchronize any sub-workflows imported by this one.

<details markdown="block">
  <summary>
    <code>manifest.json</code> format example
  </summary>
  {: .text-delta }

  <dl>
    <dt>fishers-test <code>&lt;object&gt;</code></dt>
    <dd>
      <dl>
        <dt><i>cluster_name</i> <code>&lt;str&gt;</code></dt>
        <dd>The name of the cluster-level declared in your <code>samples.csv</code> annotations for which you would like to declare a reference population for pair-wise testing.</dd>
      </dl>
    </dd>
    <dt>output-dir <code>&lt;Array&lt;Str&gt;&gt;</code></dt>
    <dd>A list representing the file-path for the location at which the workflow should save its output. If the folder does not exist, the workflow will automatically create it.</dd>
  </dl>

  ```json
  {
      "fishers-test": {
          "my_cluster": "my_population_of_interest"
      },
      "output-dir": [
          "/",
          "path",
          "to",
          "my",
          "output",
          "location"
      ]
  }
  ```
</details>
