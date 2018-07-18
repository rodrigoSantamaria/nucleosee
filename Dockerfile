FROM tiangolo/uwsgi-nginx-flask:python3.6


# Install any needed packages specified in requirements.txt
RUN pip3 install numpy
RUN pip3 install fisher 
RUN pip3 install pyBigWig

#Nginx variables
ENV STATIC_INDEX 1
#Only 1 worker to be able to share global python vars
ENV 	NGINX_WORKER_PROCESSES 1
ENV STATIC_URL /nucleosee 

#Copy my app stuff
COPY ./app /app
COPY ./uwsgi.ini /app/uwsgi.ini
COPY ./entrypoint.sh /entrypoint.sh

#RUN touch /app/genomes/tracks.txt



