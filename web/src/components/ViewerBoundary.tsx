import { Component } from "react";
import type { ReactNode } from "react";

/** A failed optional viewer chunk must not hide the catalog, results or downloads. */
export class ViewerBoundary extends Component<
  { children: ReactNode; downloadUrl: string },
  { failed: boolean }
> {
  state = { failed: false };
  static getDerivedStateFromError() {
    return { failed: true };
  }
  render() {
    if (this.state.failed)
      return (
        <div className="error" role="alert">
          The molecular viewer could not load. Reload the page to retry, or{" "}
          <a href={this.props.downloadUrl}>download the original structure</a>.
        </div>
      );
    return this.props.children;
  }
}
