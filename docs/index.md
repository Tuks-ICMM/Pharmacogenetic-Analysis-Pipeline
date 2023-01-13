---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: home
title: Welcome
permalink: /
nav_order: 1
---

<!-- START - Links, Badges and Markdown Variables -->

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg?style=for-the-badge
[snakemake-unit-tests]: https://github.com/Tuks-ICMM/Pharmacogenetic-Analysis-Pipeline/actions/workflows/snakemake-tests.yml/badge.svg

<!-- END - Links, Badges and Markdown Variables -->

![Build][snakemake-unit-tests]
![Build][docs-unit-tests]

[![CC BY 4.0][cc-by-shield]][cc-by]

<hr/>

This is the development repository for a pipeline created to perform frequency analysis on African genetic datasets.

> This work is licensed under a [Creative Commons Attribution 4.0 International License][cc-by].

# Publications and Citations

Should you wish to cite this project, please select an apropriate paper/s from the list below:

{% for post in site.categories.publication %}

> <dl>
>  <dt>Title</dt>
>  <dd><h3><a href="{{ post.permalink }}">{{ post.title }}</a></h3></dd>
>  <dt>Date</dt>
>  <dd>{{ post.date | date_to_string }}</dd>
>  <dt>DOI</dt>
>  <dd><a href="https://doi.org/10.3390/jpm12122013" target="_blank">10.3390/jpm12122013</a></dd>
> </dl>
> <hr/>
> <p>{{ post.excerpt | remove_first: "Abstract" | append: "..." }}</p>

{% endfor %}

# Roadmap:

See our [Issues tracker](https://github.com/Tuks-ICMM/Pharmacogenetic-Analysis-Pipeline/issues) on GitHub for a list of proposed features (and
known issues).

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

[![CC BY 4.0][cc-by-image]][cc-by]
