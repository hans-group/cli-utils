from rich.console import Console

from rich.theme import Theme

custom_theme = Theme({"info": "#87b8f3", "warning": "#f3bc7f", "gray": "#cdcdcd"})
console = Console(theme=custom_theme)
