export function ErrorMessage({ error }: { error: Error | null }) {
  return error ? (
    <div className="error" role="alert">
      {error.message}
    </div>
  ) : null;
}
export function Loading({
  children = "Loading…",
}: {
  children?: React.ReactNode;
}) {
  return (
    <p className="muted loading-status" role="status">
      {children}
    </p>
  );
}
export const formatCount = (value: number) => value.toLocaleString();

export function jobLabel(status: string): string {
  const labels: Record<string, string> = {
    prepared_pending_gpu_execution: "Ready to generate",
    cached_demo_completed: "Saved result",
    completed: "Complete",
    queued: "Queued",
    running: "Generating",
    dispatch_requested: "Starting",
    cancel_requested: "Cancelling",
    cancelled: "Cancelled",
    failed: "Failed",
    rejected: "Rejected",
    submission_unknown: "Submission needs review",
    submitting: "Preparing",
  };
  return labels[status] ?? status.replaceAll("_", " ");
}
