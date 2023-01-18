---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: home
title: Welcome
permalink: /
---

<details markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>

# Welcome

This is the documentation site for the Pharmacogenetics pipeline.

# Publications

<ul>
{% for post in site.categories.Publication %}
 <li>
 <span>{{ post.date | date_to_string }}</span> &nbsp; <a href="{{ post.url }}">{{ post.title }}</a>
 </li>
{% endfor %}
</ul>

# Diagrams

```mermaid

flowchart TB
input
results
subgraph val[Qaulity Control]
Validation
end

    subgraph prep[Data Preparation & Transformation]
    Collate
    Trim
    Liftover
    end

    subgraph annot[Annotation]
    Annotate
    end

    subgraph analysis[Analysis]
    Admixture
    Filter
    Transpile
    PLINK
    end

    input --> val

    %% Section links
    Validation --> prep
    prep --> annot

    Liftover --> Collate

    %% This is a fork
    annot --> Admixture
    annot --> Trim

    Trim --> Filter
    Filter --> Transpile

    Filter --> PLINK
    Transpile --> PLINK


    PLINK --> results

```
