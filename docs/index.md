---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: home
title: Welcome
permalink: /
nav_order: 1
---

<!-- START - Links, Badges and Markdown Variables -->
<div>
    <img src="https://github.com/Tuks-ICMM/Pharmacogenetic-Analysis-Pipeline/actions/workflows/snakemake-tests.yml/badge.svg"
        alt="Build status badge" />
</div>
<div>
    <a href="http://creativecommons.org/licenses/by/4.0/" target="_blank"><img
            src="https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg?style=for-the-badge"
            alt="CC BY 4.0" /></a>
    <a href="https://snakemake.readthedocs.io/" target="_blank"><img
            src="https://img.shields.io/badge/Snakemake-049474?label=Powered%20by&style=for-the-badge&logo=data:image/png;base64,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"
            alt="Powered by Snakemake" /></a>
</div>
<!-- END - Links, Badges and Markdown Variables -->
<hr />

Welcome to the <i>{{ site.title }}</i>. This is a development repository and documentation containing instruction and technical discussion on the use of {{ site.title }}. This pipeline was created by <a href="https://github.com/G-kodes" target="_blank">{% avatar G-Kodes size=15 %} Graeme Ford</a> to perform the following core functions:

- Perform VCF Spec validation
- Merge multiple VCF datasets into one psudo-dataset\*
- Perform variant frequency analysis
- Perform variant-effect-prediction

{: .normal }
> This work is licensed under a [Creative Commons Attribution 4.0 International License][cc-by].

# Publications and Citations

Should you wish to cite this project or software, please select an apropriate paper/s from the list below:

{% for post in site.categories.publication %}

<dl class="def-small">
  <dt>Title</dt>
  <dd><h3><a href="{{ post.url }}">{{ post.title }}</a></h3></dd>
  <dt>Date</dt>
  <dd>{{ post.date | date_to_string }}</dd>
  <dt>DOI</dt>
  <dd><a href="https://doi.org/10.3390/jpm12122013" target="_blank">10.3390/jpm12122013</a></dd>
  <dt>Abstract</dt>
  <dd>{{ post.content | remove_first: "Abstract" | split: site.excerpt_separator | first | strip_html }}...</dd>
</dl>

{% endfor %}


# Versioning:

We use the [SemVer](http://semver.org/) syntax to manage and maintain version numbers. For the
versions available, see the [releases on this repository
here](https://github.com/SgtPorkChops/SASDGHUB/releases).

# Acknowledgements:

Many thanks to the following individuals who have been instrumental to the success of this project.

<table>
  <tr>
    <td style="text-align:center;">
      <div>
        {% avatar G-Kodes size=100 %}
      </div>
      <h4>Author</h4>
      <hr />
      <h4>
        <strong>
          <a href="https://www.linkedin.com/in/graeme-ford/" target="_blank">
            Graeme Ford
          </a>
        </strong>
      </h4>
      <h6>
        <italic>
          <a href="https://github.com/G-kodes" target="_blank">
            G-Kodes
          </a>
        </italic>
      </h6>
    </td>
    <td style="text-align:center;">
      <div>
        {% avatar MeganHolborn size=100 %}
      </div>
      <h4>Contributor</h4>
      <hr />
      <h4>
        <strong>
          <a href="https://www.linkedin.com/in/megan-ryder-b312b0159/" target="_blank">
            Megan Holborn
          </a>
        </strong>
      </h4>
      <h6>
        <italic>
          <a href="https://github.com/orgs/Tuks-ICMM/people/MeganHolborn" target="_blank">
            MeganHolborn
          </a>
        </italic>
      </h6>
    </td>
    <td style="text-align:center;">
      <div>
        <img src="https://www.up.ac.za/media/shared/489/ZP_Images/michael-pepper-message.zp39643.jpg"
          width="100" alt="Prof. Michael S. Pepper" />
      </div>
      <h4>Supervisor</h4>
      <hr />
      <h4>
        <strong>
          <a href="https://www.up.ac.za/institute-for-cellular-and-molecular-medicine/article/2019297/professor-michael-s-pepper"
            target="_blank">Prof. Michael Pepper
          </a>
        </strong>
      </h4>
    </td>
  </tr>
  <tr>
    <td style="text-align:center;">
      <div>
          {% avatar fouriejoubert size=100 %}
      </div>
      <h4>Co-Supervisor</h4>
      <hr />
      <h4>
        <strong>
          <a href="https://www.up.ac.za/the-genomics-research-institute/article/1929131/professor-fourie-joubert"
            target="_blank">
            Prof. Fourie Joubert
          </a>
        </strong>
      </h4>
      <h6>
        <italic>
          <a href="https://github.com/orgs/Tuks-ICMM/people/fouriejoubert" target="_blank">
            fouriejoubert
          </a>
        </italic>
      </h6>
    </td>
    <td style="text-align:center;">
      <div>
        {% avatar Fatimabp size=100 %}
      </div>
      <h4>Tester</h4>
      <hr />
      <h4>
        <strong>
          <a href="https://www.linkedin.com/in/fatima-barmania-a1201238/" target="_blank">
            Fatima Barmania-Faruk
          </a>
        </strong>
      </h4>
      <h6>
        <italic>
          <a href="https://github.com/orgs/Tuks-ICMM/people/Fatimabp" target="_blank">
            Fatimabp
          </a>
        </italic>
      </h6>
    </td>
    <td style="text-align:center;">
      <div>
        {% avatar sarahsaraht size=100 %}
      </div>
      <h4>Tester</h4>
      <hr />
      <h4>
        <strong>
          <a href="https://upload.wikimedia.org/wikipedia/commons/c/cd/Portrait_Placeholder_Square.png" target="_blank">
            Sarah Turner
          </a>
        </strong>
      </h4>
      <h6>
        <italic>
          <a href="https://github.com/orgs/Tuks-ICMM/people/sarahsaraht" target="_blank">
            sarahsaraht
          </a>
        </italic>
      </h6>
    </td>
  </tr>
</table>

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png

[![CC BY 4.0][cc-by-image]][cc-by]
