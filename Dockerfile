FROM python:3.8-slim

WORKDIR /usr/src/app

COPY . /usr/src/app

RUN pip3 install matplotlib sympy numpy pexpect

CMD ["python3", "discrete_logarithm.py"] 