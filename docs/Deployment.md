# Dashboard Deployment

## Systemd Service

To keep the Plotly Dash app running in the background and restart it automatically, install the systemd unit provided in the repository:

```bash
sudo cp dashboard/gvmagdb-dashboard.service /etc/systemd/system/
sudo systemctl daemon-reload
sudo systemctl enable --now gvmagdb-dashboard.service
```

The service assumes the repository lives at `/home/fschulz/dev/gvmagDB` and runs with the Pixi-managed interpreter (`.pixi/envs/default/bin/python`). Adjust `WorkingDirectory`, `ExecStart`, or `User` if your paths differ.

Key commands:

- `sudo systemctl status gvmagdb-dashboard.service`
- `sudo journalctl -u gvmagdb-dashboard.service -f`
- `sudo systemctl restart gvmagdb-dashboard.service`

## Tailscale Funnel

Expose the dashboard publicly (over HTTPS) via Tailscale Funnel:

```bash
sudo tailscale funnel --bg --https=443 --set-path=/ 127.0.0.1:8050/
```

Verify:

```bash
tailscale funnel status
```

To disable: `sudo tailscale funnel --https=443 off`

## Analytics Cache Refresh

The Dash app reads cached analytics tables from `artifacts/analytics/` (see `pixi run dashboard-cache`). Refresh the cache after ingestion or whenever the dataset updates:

```bash
pixi run dashboard-cache
sudo systemctl restart gvmagdb-dashboard.service
```

Automate via cron or a scheduled systemd timer if needed.
