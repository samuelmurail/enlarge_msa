"""Console script for enlarge_msa."""

import enlarge_msa

import typer
from rich.console import Console

app = typer.Typer()
console = Console()


@app.command()
def main():
    """Console script for enlarge_msa."""
    console.print(
        "Replace this message by putting your code into " "enlarge_msa.cli.main"
    )
    console.print("See Typer documentation at https://typer.tiangolo.com/")


if __name__ == "__main__":
    app()
