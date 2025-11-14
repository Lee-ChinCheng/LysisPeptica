# 1. Start from a base image with Conda
FROM continuumio/miniconda3:latest

# 2. Set the working directory inside the container
WORKDIR /app

# 3. Copy files needed to build the environment
COPY environment.yml .

# 4. Create the Conda environment
RUN conda env create -f environment.yml

# 5. Copy the models, encoding function, and prediction script
COPY models/ /app/models/
COPY config.py .
COPY func.py .
COPY model.py .
COPY predictor.py .
COPY tokenizer.json .
COPY UniRef50_len_0_50loss.csv .

# 6. Create directories for I/O
RUN mkdir /app/data
RUN mkdir /app/output

# 7. Set the default command (ENTRYPOINT)
# This now passes the standard I/O paths.
ENTRYPOINT ["conda", "run", "-n", "py310", "python", "predictor.py", \
            "--input", "/app/data/input.fasta", \
            "--output", "/app/output/predictions.csv"]
