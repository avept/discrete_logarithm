FROM python:3.11-slim

WORKDIR /script
COPY script/. .

RUN pip install --upgrade pip
RUN pip install numpy
RUN pip install sympy
RUN pip install argparse

CMD ["python", "discrete_logarithm.py"]
