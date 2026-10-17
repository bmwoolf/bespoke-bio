from fastapi import FastAPI, HTTPException

from bespoke_bio.models import TranslateRequest, TranslateResponse
from bespoke_bio.pipeline import translate_case
from bespoke_bio.translation import VariantApplicationError

app = FastAPI(
    title="Bespoke Bio",
    summary="Translate expressed rare variants into target proteins",
    version="0.1.0",
)


@app.get("/healthz")
def healthcheck() -> dict[str, str]:
    return {"status": "ok"}


@app.post("/v1/translate", response_model=TranslateResponse)
def translate(request: TranslateRequest) -> TranslateResponse:
    try:
        return translate_case(request)
    except VariantApplicationError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error
