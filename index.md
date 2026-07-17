---
title: Home
layout: default
nav_order: 0
---

![Banner](./docs/imgs/banner_sasmoca_light.png)

# SAS MoCa Documentation

🚧 Work in progress

## Navigation

<ul>
{% assign sorted_pages = site.pages | sort: "nav_order" %}
{% for page in sorted_pages %}
  {% if page.layout == "default" or page.layout == "home" %}
    <li><a href="{{ page.url }}">{{ page.title }}</a></li>
  {% endif %}
{% endfor %}
</ul>

## Documentation Pages

| Page                                       | Description                  |
| ------------------------------------------ | ---------------------------- |
| [Getting Started](docs/getting_started.md) | Installation and first steps |
