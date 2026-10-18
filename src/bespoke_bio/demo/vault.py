"""Deployment helper: fetch one credential using the VM identity, without logging its value."""

import argparse
import os
from pathlib import Path
from tempfile import NamedTemporaryFile

from .client import DesignService


def materialize(value, destination):
    """Validate before atomically replacing a private runtime credential file."""
    destination = Path(destination)
    with NamedTemporaryFile(mode="w", dir=destination.parent, delete=False) as stream:
        temporary = Path(stream.name)
        try:
            os.chmod(temporary, 0o600)
            stream.write(value)
            stream.flush()
            os.fsync(stream.fileno())
            DesignService.from_config(temporary).check_available()
            temporary.replace(destination)
        finally:
            temporary.unlink(missing_ok=True)


def main():
    from azure.identity import ManagedIdentityCredential
    from azure.keyvault.secrets import SecretClient

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vault-url", required=True)
    parser.add_argument("--secret", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    with (
        ManagedIdentityCredential() as credential,
        SecretClient(vault_url=args.vault_url, credential=credential) as client,
    ):
        value = client.get_secret(args.secret).value
        materialize(value, args.output)


if __name__ == "__main__":
    main()
