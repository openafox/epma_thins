---
layout: default
title: Home
---

# EPMA Thin Film Analysis
{: style="text-align: center;"}
A Python package still under development based on [GMRFilm](https://github.com/openafox/gmrfilm) by R.A. Waldo.

Some of the code is a direct port from the original fortran code but much is not and is refactored to a OOP style.

The working code can be found on [GitHub](https://github.com/openafox/epma_thins).

If you would like to contribute please see [contributing](/contributing). I really appreciate any help making this a reality.

Some of the math has been checked from the initial publications (see the Math category).

  <div class="tags-expo-list" style="text-align: center;">
    {% for tag in site.categories %}
    <a href="/blog/categories#{{ tag[0] | slugify }}" class="post-tag">{{ tag[0] }}</a>
    {% endfor %}
  </div>
