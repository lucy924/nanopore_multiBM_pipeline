## Basic output
**Sample Name:** {{repo.name}}<br>
**URL:** {{repo.url}}

| Category | Result |
|---|---|
{% for entry in repo.results %}
|{{entry[0]}}|{{entry[1]}}|
{% endfor %}

