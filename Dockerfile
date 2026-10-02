FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    HOST=0.0.0.0 \
    PORT=8080

WORKDIR /srv

COPY app ./app
COPY tests ./tests
COPY scripts ./scripts

# No third-party runtime dependencies; this also fails the build if any
# source file has a syntax error.
RUN python -m compileall -q app scripts

EXPOSE 8080

HEALTHCHECK --interval=5s --timeout=3s --start-period=3s --retries=10 \
    CMD python -c "import json,urllib.request,sys; \
r=urllib.request.urlopen('http://127.0.0.1:8080/healthz',timeout=3); \
sys.exit(0 if r.status==200 and json.load(r)['status']=='ok' else 1)"

CMD ["python", "-m", "app.main"]
