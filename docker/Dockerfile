FROM python:3.12
WORKDIR /usr/local

COPY requirements.txt .
RUN apt-get install git
RUN pip install --no-cache-dir -r requirements.txt
RUN pip install git+https://github.com/cellgeni/cloupe.git
RUN git clone https://github.com/hassansaei/snapseed.git
RUN pip install snapseed/

CMD ["/bin/bash"]
