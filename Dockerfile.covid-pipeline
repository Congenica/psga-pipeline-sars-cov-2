FROM python:3.8

# pin the versions of pip and poetry
ENV PIP_VERSION=20.2.4 \
    POETRY_VERSION=1.1.4

RUN apt-get update && \
    apt-get install --yes --no-install-recommends build-essential procps libpq-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY . .

# used by our reports
RUN wget -P /app/data/pango-designation \
    https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineages.csv \
    https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt

RUN pip install --progress-bar=off --no-cache-dir --upgrade pip==${PIP_VERSION} && \
    pip install --progress-bar=off --no-cache-dir poetry==${POETRY_VERSION} && \
    # make poetry install everything system wide, no need for a virtualenv in docker
    poetry config virtualenvs.create false && \
    poetry install --no-interaction --no-ansi --no-dev
