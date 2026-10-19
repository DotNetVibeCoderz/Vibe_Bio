namespace Biotech.Net.Workflows;

/// <summary>Something that happened during a run; delivered through <see cref="IProgress{T}"/> or <see cref="Workflow.RunStreamingAsync"/>.</summary>
/// <param name="Timestamp">When it happened, from the run's <see cref="TimeProvider"/>.</param>
public abstract record WorkflowEvent(DateTimeOffset Timestamp);

/// <summary>The run started.</summary>
public sealed record RunStartedEvent(DateTimeOffset Timestamp, string RunId, int JobCount) : WorkflowEvent(Timestamp);

/// <summary>A job was up to date (or is a pure aggregate) and did not run.</summary>
public sealed record JobSkippedEvent(DateTimeOffset Timestamp, string JobId, string Reason) : WorkflowEvent(Timestamp);

/// <summary>An attempt of a job started.</summary>
public sealed record JobStartedEvent(DateTimeOffset Timestamp, string JobId, int Attempt, string Reason) : WorkflowEvent(Timestamp);

/// <summary>A log line from a job (its own messages, or a line of a tool's stdout/stderr).</summary>
public sealed record JobLogEvent(DateTimeOffset Timestamp, string JobId, string Message) : WorkflowEvent(Timestamp);

/// <summary>An attempt failed and the job will be retried after <paramref name="Delay"/>.</summary>
public sealed record JobRetryingEvent(DateTimeOffset Timestamp, string JobId, int FailedAttempt, TimeSpan Delay, string Error) : WorkflowEvent(Timestamp);

/// <summary>A job finished and its outputs were moved into place.</summary>
public sealed record JobSucceededEvent(DateTimeOffset Timestamp, string JobId, TimeSpan Duration) : WorkflowEvent(Timestamp);

/// <summary>A job failed for good (all attempts used, or not retryable).</summary>
public sealed record JobFailedEvent(DateTimeOffset Timestamp, string JobId, string Error) : WorkflowEvent(Timestamp);

/// <summary>A job was cancelled while running.</summary>
public sealed record JobCancelledEvent(DateTimeOffset Timestamp, string JobId) : WorkflowEvent(Timestamp);

/// <summary>Jobs finished so far (run, skipped, failed or cancelled) out of <paramref name="Total"/>.</summary>
public sealed record RunProgressEvent(DateTimeOffset Timestamp, int Finished, int Total) : WorkflowEvent(Timestamp);

/// <summary>The run ended; the report is also written to the state directory.</summary>
public sealed record RunCompletedEvent(DateTimeOffset Timestamp, WorkflowRunReport Report) : WorkflowEvent(Timestamp);
