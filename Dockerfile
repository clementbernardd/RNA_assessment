FROM python:3.8
WORKDIR app
COPY DEPENDENCY /app/
RUN pip install -r DEPENDENCY
COPY . .
CMD ["python", "-m", "example"]
